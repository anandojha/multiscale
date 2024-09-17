from openmm.app import StateDataReporter
from openmm.app import DCDReporter
from openmm.app import PDBReporter
from openmm.app import Simulation
import openmm.openmm as openmm
from openmm.app import PDBFile
from openmm import unit
from openmm import app
import openmm as mm
import sys
import os

from openmm.app import StateDataReporter, DCDReporter, Simulation, PDBFile
import openmm as mm
from openmm import unit
import sys
import os

def minimization_simulation(pdb_file, system_xml, output_prefix, platform_name='CUDA', temperature=300, friction_coeff=1, time_step=0.002, max_minimize_steps=100000, md_steps=100000, report_interval=10000):
    """
    Run an OpenMM molecular dynamics simulation from a given PDB file and serialized system XML.

    Parameters:
    pdb_file (str): Path to the PDB file for setting up initial positions.
    system_xml (str): Path to the serialized system XML file.
    output_prefix (str): Prefix for output files.
    platform_name (str): Platform to use for the simulation ('CUDA', 'OpenCL', 'CPU').
    temperature (float): Simulation temperature in Kelvin.
    friction_coeff (float): Friction coefficient for Langevin integrator in ps^-1.
    time_step (float): Time step for the integrator in ps.
    max_minimize_steps (int): Maximum steps for energy minimization.
    md_steps (int): Number of steps for MD simulation.
    report_interval (int): Interval at which simulation data is reported.

    Returns:
    None
    """
    # Load the PDB file and serialized system XML
    pdb = PDBFile(pdb_file)
    with open(system_xml, 'r') as f:
        system = mm.XmlSerializer.deserialize(f.read())
    # Set up the integrator for the simulation
    integrator = mm.LangevinIntegrator(temperature * unit.kelvin, friction_coeff / unit.picosecond, time_step * unit.picoseconds)
    # Choose the platform for simulation
    platform = mm.Platform.getPlatformByName(platform_name)
    # Set up the simulation object
    simulation = Simulation(pdb.topology, system, integrator, platform)
    # Set initial positions from the PDB file
    simulation.context.setPositions(pdb.positions)
    # Minimize the energy of the system
    print(f"Minimizing energy for {max_minimize_steps} steps...")
    simulation.minimizeEnergy(maxIterations=max_minimize_steps)
    # Print the minimized energy
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    print(f"Minimized energy: {state.getPotentialEnergy()}")
    # Save the minimized structure as a serialized XML file
    minimized_system_xml = mm.XmlSerializer.serialize(system)
    with open(f'{output_prefix}_minimized_system.xml', 'w') as f:
        f.write(minimized_system_xml)
    print(f"Minimized system saved as '{output_prefix}_minimized_system.xml'.")
    # Save the minimized structure as a PDB file
    with open(f'{output_prefix}_minimized_structure.pdb', 'w') as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)
    print(f"Minimized structure saved as '{output_prefix}_minimized_structure.pdb'.")
    # Add reporters to save simulation data
    simulation.reporters.append(StateDataReporter(sys.stdout, report_interval, step=True, potentialEnergy=True, temperature=True))
    simulation.reporters.append(DCDReporter(f'{output_prefix}_trajectory.dcd', report_interval))
    # Run the MD simulation
    print("Running MD simulation...")
    simulation.step(md_steps)
    # Save the final frame of the MD simulation as a PDB file
    final_state = simulation.context.getState(getPositions=True)
    #with open(f'{output_prefix}_final_structure.pdb', 'w') as f:
        #PDBFile.writeFile(simulation.topology, final_state.getPositions(), f)
    #print(f"Final structure saved as '{output_prefix}_final_structure.pdb'.")

def run_minimization_simulation(output_prefix='output', platform_name='CUDA'):
    # List all files in the current directory
    files = os.listdir('.')
    # Initialize variables to store the identified file names
    pdb_file = None
    system_xml = None
    # Identify the .pdb and .xml files
    for file in files:
        if file.endswith('.pdb'):
            pdb_file = file
        elif file.endswith('.xml'):
            system_xml = file
    # Check if both files are identified
    if pdb_file and system_xml:
        print(f"Identified PDB file: {pdb_file}")
        print(f"Identified XML file: {system_xml}")
        # Call the minimization_simulation function
        minimization_simulation(
            pdb_file=pdb_file,
            system_xml=system_xml,
            output_prefix=output_prefix,
            platform_name=platform_name
        )
    else:
        raise FileNotFoundError("Required .pdb or .xml files not found in the current directory.")

def cleanup_rename_files():
    # List all files in the current directory
    files = os.listdir('.')
    # Initialize variables to store the identified original file names
    original_pdb_file = None
    original_xml_file = None
    generated_files = []
    # Identify the original .pdb and .xml files
    for file in files:
        if file.endswith('.pdb') and 'output' not in file:
            original_pdb_file = file
        elif file.endswith('.xml') and 'output' not in file:
            original_xml_file = file
    # Identify generated files (those starting with 'output')
    for file in files:
        if file.startswith('output') and (file.endswith('.xml') or file.endswith('.dcd')):
            generated_files.append(file)
    # Ensure necessary files were identified
    if not original_pdb_file or not original_xml_file:
        raise FileNotFoundError("Original .pdb or .xml file not found in the current directory.")
    # Delete unnecessary generated files
    for gen_file in generated_files:
        os.remove(gen_file)
        print(f"Deleted generated file: {gen_file}")
    # Rename the original .pdb file
    original_pdb_base_name = original_pdb_file[:-4]  # Remove '.pdb' extension
    new_original_pdb_name = f"{original_pdb_base_name}_from_metaD_SMD.pdb"
    os.rename(original_pdb_file, new_original_pdb_name)
    print(f"Renamed original PDB file from {original_pdb_file} to {new_original_pdb_name}")
    # Rename the generated file output_minimized_structure.pdb to the original .pdb filename
    generated_minimized_pdb = 'output_minimized_structure.pdb'
    new_generated_pdb_name = original_pdb_file
    os.rename(generated_minimized_pdb, new_generated_pdb_name)
    print(f"Renamed {generated_minimized_pdb} to {new_generated_pdb_name}")


run_minimization_simulation()
cleanup_rename_files()


























