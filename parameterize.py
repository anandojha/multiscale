"""
Create complex system with ESPALOMA.
 - ESPALOMA is used as the default force field to parameterize both small molecules and proteins.
 - ESPALOMA and AMBER FORCE FIELD are used to first parameterize the ligand and protein, respectively, to create the complex with solvation and ions.
 - The protein and the ligand from the solvated complex is then extracted and re-parameterized with ESPALOMA.
 - Note that the script works when the oe_license.txt is exported 
export OE_LICENSE="/home/aaojha/oe_license.txt"

To install ESPALOMA:
conda create -n ESPALOMA python=3.10 --y
conda activate ESPALOMA
conda install conda-forge::openff-toolkit --y
conda install conda-forge::openmmforcefields --y
conda install -c conda-forge "espaloma=0.3.2" --y
conda install openeye::openeye-toolkits --y
"""

# Import necessary modules from OpenMM, OpenFF, RDKit, tqdm, mdtraj, numpy, warnings, requests, zipfile, logging, shutil, and os
from openmm import MonteCarloBarostat, LangevinMiddleIntegrator, XmlSerializer
from openmmforcefields.generators import SystemGenerator
from openff.toolkit.topology import Molecule
import openmm.unit as unit
import openmm.app as app
from rdkit import Chem
from tqdm import tqdm
import mdtraj as md
import numpy as np
import warnings
import requests
import zipfile
import logging
import shutil
import os

# Suppress all warnings to keep output clean
warnings.filterwarnings("ignore")

# Set up logging to record information in a file named 'parameterize.log'
logging.basicConfig(filename='parameterize.log', level=logging.INFO)
_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)

def get_forcefields_dir():
    """
    Download and extract forcefields.

    - Downloads forcefields.zip from the provided URL.
    - Extracts the contents of the zip file into the current directory.
    - Removes any existing 'forcefields' directory before extraction to avoid conflicts.
    """
    current_directory = os.getcwd()  # Get current working directory
    forcefields_directory = os.path.join(current_directory, "forcefields")  # Define forcefields directory path
    if os.path.exists(forcefields_directory):
        shutil.rmtree(forcefields_directory)  # Remove existing directory
        print(f"Removed existing directory: {forcefields_directory}")
    url = "https://zenodo.org/record/12797986/files/forcefields.zip?download=1"  # URL to download forcefields
    zip_path = os.path.join(current_directory, "forcefields.zip")  # Define path for the downloaded zip file
    response = requests.get(url, stream=True)  # Send GET request to download the file
    total_size = int(response.headers.get('content-length', 0))  # Get total size of the file for progress bar
    with open(zip_path, "wb") as file, tqdm(desc="Downloading forcefields.zip", total=total_size, unit='B', unit_scale=True, unit_divisor=1024) as bar:
        for data in response.iter_content(chunk_size=1024):
            file.write(data)  # Write data to file in chunks
            bar.update(len(data))  # Update progress bar
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(current_directory)  # Extract all files to the current directory
    nested_forcefields_directory = os.path.join(forcefields_directory, "forcefields")
    if os.path.exists(nested_forcefields_directory):
        for item in os.listdir(nested_forcefields_directory):
            shutil.move(os.path.join(nested_forcefields_directory, item), forcefields_directory)  # Move nested items to main directory
        shutil.rmtree(nested_forcefields_directory)  # Remove nested directory
    os.remove(zip_path)  # Remove the zip file after extraction

def create_complex(protein_file, ligand_file, forcefield_files, water_model, solvent_padding, pdb_template, ionic_strength, hmass, temperature, pressure, pme_tol, barostat_period, system_forcefield, friction, stepsize, maxiterations, serialized_xml, state_xml, minimized_pdb, integrator_xml, processed_pdb):
    """
    Create a complex of protein and ligand, solvate, and parameterize with ESPALOMA.

    Parameters:
    - protein_file: Path to the protein PDB file.
    - ligand_file: Path to the ligand SDF file.
    - forcefield_files: List of forcefield XML files.
    - water_model: Water model to use (e.g., 'tip3p').
    - solvent_padding: Padding for solvation in Angstroms.
    - pdb_template: PDB file where pritein and ligand atoms are saved. 
    - ionic_strength: Ionic strength for solvation.
    - hmass: Hydrogen mass in atomic mass units (amu).
    - temperature: Temperature in Kelvin.
    - pressure: Pressure in atmosphere.
    - pme_tol: Ewald error tolerance for PME.
    - barostat_period: Frequency of barostat updates.
    - system_forcefield: Force field for the system.
    - friction: Friction coefficient for Langevin integrator.
    - stepsize: Integration step size in femtoseconds.
    - maxiterations: Maximum iterations for energy minimization.
    - serialized_xml: Path to save serialized system XML.
    - state_xml: Path to save serialized state XML.
    - minimized_pdb: Path to save minimized PDB file.
    - integrator_xml: Path to save serialized integrator XML.
    - processed_pdb: Processed PDB file for MD simulation.
    """
    # Load ligand from SDF file
    ext = os.path.splitext(ligand_file)[-1].lower()  # Get file extension
    assert ext == '.sdf', f'Ligand file format must be SDF but got {ext}'  # Ensure file is in SDF format
    suppl = Chem.SDMolSupplier(ligand_file)  # Load molecules from SDF
    mols = [x for x in suppl if x is not None]  # Filter out None molecules
    mol = mols[0]  # Select the first molecule
    mol.SetProp("_Name", "MOL")  # Set molecule name (optional)
    offmol = Molecule.from_rdkit(mol)  # Convert RDKit molecule to OpenFF molecule
    ligand_positions = offmol.conformers[0] / 10  # Get ligand positions and convert to nanometers

    # Load protein from PDB file
    with open(protein_file, 'r') as f:
        protein = app.PDBFile(f)  # Load protein PDB file

    # Merge topologies of protein and ligand
    protein_topology = protein.topology  # Get protein topology
    protein_positions = protein.positions  # Get protein positions
    ligand_topology = offmol.to_topology().to_openmm()  # Convert ligand topology to OpenMM format
    protein_md_topology = md.Topology.from_openmm(protein_topology)  # Convert protein topology to MDTraj format
    ligand_md_topology = md.Topology.from_openmm(ligand_topology)  # Convert ligand topology to MDTraj format
    complex_md_topology = protein_md_topology.join(ligand_md_topology)  # Join protein and ligand topologies
    complex_topology = complex_md_topology.to_openmm()  # Convert merged topology back to OpenMM format

    # Ensure the number of atoms is consistent after merging
    n_atoms_total = complex_md_topology.n_atoms
    n_atoms_protein = protein_md_topology.n_atoms
    n_atoms_ligand = ligand_md_topology.n_atoms
    assert n_atoms_total == n_atoms_protein + n_atoms_ligand, "Mismatch in atom numbers after merging."

    # Combine positions of protein and ligand
    complex_positions = unit.Quantity(np.zeros([n_atoms_total, 3]), unit=unit.nanometers)  # Initialize complex positions
    protein_positions_nm = unit.Quantity(np.array(protein.positions / unit.nanometers), unit.nanometers)  # Convert protein positions to nanometers
    complex_positions[:n_atoms_protein, :] = protein_positions_nm  # Assign protein positions
    ligand_positions_nm = unit.Quantity(ligand_positions, unit.nanometers)  # Convert ligand positions to nanometers
    complex_positions[n_atoms_protein:n_atoms_protein + n_atoms_ligand, :] = ligand_positions_nm  # Assign ligand positions

    # Initialize system generator with force field parameters
    forcefield_kwargs = {'removeCMMotion': True, 'ewaldErrorTolerance': pme_tol, 'constraints': app.HBonds, 'rigidWater': True, 'hydrogenMass': hmass * unit.amu}
    periodic_forcefield_kwargs = {'nonbondedMethod': app.PME}
    barostat = MonteCarloBarostat(pressure * unit.atmosphere, temperature * unit.kelvin, barostat_period)  # Initialize barostat
    system_generator = SystemGenerator(forcefields=forcefield_files, forcefield_kwargs=forcefield_kwargs, periodic_forcefield_kwargs=periodic_forcefield_kwargs, barostat=barostat, small_molecule_forcefield=system_forcefield, molecules=offmol, cache=None)

    # Solvate the system
    modeller = app.Modeller(complex_topology, complex_positions)  # Initialize modeller with complex topology and positions
    modeller.addSolvent(system_generator.forcefield, model=water_model, padding=solvent_padding * unit.angstroms, ionicStrength=ionic_strength * unit.molar)  # Add solvent
    solvated_topology = modeller.getTopology()  # Get solvated topology
    solvated_positions = modeller.getPositions()  # Get solvated positions
    solvated_system = system_generator.create_system(solvated_topology)  # Create system with solvated topology

    # Export the complex to a PDB file
    app.PDBFile.writeFile(solvated_topology, solvated_positions, file=open(pdb_template, 'w'))

    # Regenerate the system with ESPALOMA
    regenerate_espaloma_system(system_generator, solvated_topology, solvated_positions, temperature=temperature, friction=friction, stepsize=stepsize, maxiterations=maxiterations, serialized_xml=serialized_xml, state_xml=state_xml, minimized_pdb=minimized_pdb, integrator_xml=integrator_xml)

    # Process the PDB file for simulation 
    process_pdb_file(pdb_template=pdb_template, pdb_to_process=minimized_pdb, processed_pdb=processed_pdb)

def regenerate_espaloma_system(system_generator, solvated_topology, solvated_positions, temperature, friction, stepsize, maxiterations, serialized_xml, state_xml, minimized_pdb, integrator_xml):
    """
    Regenerate the system with ESPALOMA.

    Parameters:
    - system_generator: The system generator.
    - solvated_topology: Solvated topology of the complex.
    - solvated_positions: Solvated positions of the complex.
    - temperature: Temperature in Kelvin.
    - friction: Friction coefficient for Langevin integrator.
    - stepsize: Integration step size in femtoseconds.
    - maxiterations: Maximum iterations for energy minimization.
    - serialized_xml: Path to save serialized system XML.
    - state_xml: Path to save serialized state XML.
    - minimized_pdb: Path to save minimized PDB file.
    - integrator_xml: Path to save serialized integrator XML.
    """
    # Convert solvated topology to MDTraj format and identify protein chains
    mdtop = md.Topology.from_openmm(solvated_topology)
    chain_indices = [chain.index for chain in solvated_topology.chains()]
    protein_chain_indices = [chain_index for chain_index in chain_indices if mdtop.select(f"protein and chainid == {chain_index}").any()]

    # Create new topology and copy chains
    new_solvated_topology = app.Topology()
    new_solvated_topology.setPeriodicBoxVectors(solvated_topology.getPeriodicBoxVectors())
    new_atoms = {}
    chain_counter = 0
    for chain in solvated_topology.chains():
        new_chain = new_solvated_topology.addChain(chain.id)
        if chain.index in protein_chain_indices:
            resname = f'XX{chain_counter:01d}'  # Assign unique residue name for each protein chain
            resid = '1'
            chain_counter += 1
            new_residue = new_solvated_topology.addResidue(resname, new_chain, resid)
        for residue in chain.residues():
            if residue.chain.index not in protein_chain_indices:
                new_residue = new_solvated_topology.addResidue(residue.name, new_chain, residue.id)
            for atom in residue.atoms():
                new_atom = new_solvated_topology.addAtom(atom.name, atom.element, new_residue, atom.id)
                new_atoms[atom] = new_atom
    for bond in solvated_topology.bonds():
        if bond[0] in new_atoms and bond[1] in new_atoms:
            new_solvated_topology.addBond(new_atoms[bond[0]], new_atoms[bond[1]])

    # Save the complex with ESPALOMA parameterization to PDB file
    complex_espaloma_filename = f"complex_espaloma.pdb"
    app.PDBFile.writeFile(new_solvated_topology, solvated_positions, file=open(complex_espaloma_filename, 'w'))

    # Split protein chains into separate PDB files
    protein_espaloma_filenames = []
    for chain_index in protein_chain_indices:
        t = md.load_pdb(complex_espaloma_filename)
        indices = t.topology.select(f"chainid == {chain_index}")
        t.atom_slice(indices).save_pdb(f"complex_espaloma-{chain_index}.pdb")
        protein_espaloma_filenames.append(f"complex_espaloma-{chain_index}.pdb")

    # Load protein molecules and add to system generator template
    protein_molecules = [Molecule.from_file(protein_filename) for protein_filename in protein_espaloma_filenames]
    system_generator.template_generator.add_molecules(protein_molecules)

    # Create new solvated system and minimize energy
    new_solvated_system = system_generator.create_system(new_solvated_topology)
    minimize_system(new_solvated_topology, solvated_positions, new_solvated_system, temperature=temperature, friction=friction, stepsize=stepsize, maxiterations=maxiterations, serialized_xml=serialized_xml, state_xml=state_xml, minimized_pdb=minimized_pdb, integrator_xml=integrator_xml)

def minimize_system(topology, positions, system, temperature, friction, stepsize, maxiterations, serialized_xml, state_xml, minimized_pdb, integrator_xml):
    """
    Minimize the energy of the solvated system.

    Parameters:
    - topology: System topology.
    - positions: Atom positions.
    - system: OpenMM system.
    - temperature: Temperature in Kelvin.
    - friction: Friction coefficient for Langevin integrator.
    - stepsize: Integration step size in femtoseconds.
    - maxiterations: Maximum iterations for energy minimization.
    - serialized_xml: Path to save serialized system XML.
    - state_xml: Path to save serialized state XML.
    - minimized_pdb: Path to save minimized PDB file.
    - integrator_xml: Path to save serialized integrator XML.
    """
    integrator = LangevinMiddleIntegrator(temperature * unit.kelvin, friction/unit.picosecond, stepsize * unit.femtoseconds)  # Initialize Langevin integrator
    simulation = app.Simulation(topology, system, integrator)  # Initialize simulation with system and integrator
    simulation.context.setPositions(positions)  # Set atom positions
    simulation.minimizeEnergy(maxIterations=maxiterations)  # Minimize energy
    export_system(system, simulation, serialized_xml=serialized_xml, state_xml=state_xml, minimized_pdb=minimized_pdb, integrator_xml=integrator_xml)

def export_system(system, simulation, serialized_xml, state_xml, minimized_pdb, integrator_xml):
    """
    Serialize and export the system and its state.

    Parameters:
    - system: OpenMM system.
    - simulation: OpenMM simulation object.
    - serialized_xml: Path to save serialized system XML.
    - state_xml: Path to save serialized state XML.
    - minimized_pdb: Path to save minimized PDB file.
    - integrator_xml: Path to save serialized integrator XML.
    """
    state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True)  # Get simulation state
    with open(serialized_xml, "w") as wf:
        xml = XmlSerializer.serialize(system)  # Serialize system
        wf.write(xml)
    with open(state_xml, "w") as wf:
        xml = XmlSerializer.serialize(state)  # Serialize state
        wf.write(xml)
    with open(minimized_pdb, "w") as wf:
        app.PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), file=wf, keepIds=True)  # Save minimized positions to PDB
    with open(integrator_xml, "w") as wf:
        xml = XmlSerializer.serialize(simulation.integrator)  # Serialize integrator
        wf.write(xml)

def process_pdb_file(pdb_template, pdb_to_process, processed_pdb):
    """
    Combine columns from two PDB files and output the result.

    This function reads two PDB files, combines the first 27 columns from the
    template PDB file with the remaining columns from the PDB file to process,
    ensures the atom names and counts match, and writes the combined data to a
    new output PDB file. The crystal information and headers from the template
    PDB file are preserved.

    Parameters:
    - pdb_template: str
        Path to the template PDB file which provides the first 27 columns of the ATOM/HETATM lines.
    - pdb_to_process: str
        Path to the PDB file that provides the columns after the 27th column of the ATOM/HETATM lines.
    - processed_pdb: str
        Path to save the combined output PDB file.

    Raises:
    - ValueError: If the atom names or the number of atoms do not match between the two PDB files.
    """
    
    pdb_template_lines = []
    pdb_to_process_lines = []
    
    # Read the lines of the PDB files
    with open(pdb_template, 'r') as f1:
        pdb_template_lines = f1.readlines()
    with open(pdb_to_process, 'r') as f2:
        pdb_to_process_lines = f2.readlines()
    
    combined_lines = []
    header_lines = []
    
    # Extract header and crystal information from pdb_template
    for line in pdb_template_lines:
        if line.startswith(('ATOM', 'HETATM')):
            break
        header_lines.append(line)
    
    for line1, line2 in zip(pdb_template_lines, pdb_to_process_lines):
        if line1.startswith(('ATOM', 'HETATM')) and line2.startswith(('ATOM', 'HETATM')):
            # Check that the atom names and number of atoms match
            atom_name1 = line1[12:16].strip()
            atom_name2 = line2[12:16].strip()
            if atom_name1 != atom_name2:
                raise ValueError(f"Atom names do not match: {atom_name1} vs {atom_name2}")
            
            # Combine columns 1-27 from line1 with the rest from line2
            combined_line = line1[:27] + line2[27:]
            combined_lines.append(combined_line)
    
    # Write the combined lines to the output PDB file
    with open(processed_pdb, 'w') as out_f:
        # Write header lines first
        out_f.writelines(header_lines)
        # Write combined ATOM/HETATM lines
        out_f.writelines(combined_lines)


# Download the forcefields directory
#get_forcefields_dir()

# Create complex system with given parameters
create_complex(protein_file='protein.pdb', ligand_file='ligand.sdf', forcefield_files=['amber/ff14SB.xml','amber/tip3p_standard.xml','amber/tip3p_HFE_multivalent.xml'], water_model='tip3p', solvent_padding=9.0, pdb_template='complex.pdb', ionic_strength=0.15, hmass=1.0, temperature=300.0, pressure=1.0, pme_tol=2.5e-04, barostat_period=50, system_forcefield='espaloma-0.3.2.pt', friction=1.0, stepsize=2, maxiterations=10000000, serialized_xml='complex_serialized.xml', state_xml='complex_state.xml', minimized_pdb='complex_minimized.pdb', integrator_xml='complex_integrator.xml', processed_pdb='complex_SEEKR.pdb')
