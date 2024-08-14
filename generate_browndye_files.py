from parmed import openmm as parmed_openmm
from simtk.openmm import XmlSerializer
from simtk.openmm.app import PDBFile
import os

def get_pqr(xml_file, pdb_file, output_file):
    """
    Load system from an XML file, coordinates from a PDB file, filter out unwanted residues,
    and save the filtered structure as a PQR file.
    
    Args:
        xml_file (str): Path to the XML file containing system data.
        pdb_file (str): Path to the PDB file containing coordinates.
        output_file (str): Path where the output PQR file will be saved.
    """
    residues = ["HOH", "WAT", "NA", "CL"]
    # Load the system from the XML file
    with open(xml_file, 'r') as file:
        system = XmlSerializer.deserialize(file.read())
    # Load the coordinates from the PDB file
    pdb = PDBFile(pdb_file)
    positions = pdb.getPositions(asNumpy=True)
    # Load the structure using ParmEd
    structure = parmed_openmm.load_topology(pdb.topology, system)
    structure.positions = positions
    # Function to filter out unwanted residues from the structure
    def filter_structure(structure, residues):
        filtered_indices = []
        for residue in structure.residues:
            if residue.name not in residues:
                filtered_indices.extend(atom.idx for atom in residue.atoms)
        # Create a new structure with filtered atoms
        filtered_structure = structure[filtered_indices]
        return filtered_structure
    # Filter out unwanted residues
    filtered_structure = filter_structure(structure, residues)
    # Save the filtered structure as a PQR file
    if os.path.exists(output_file):
        os.remove(output_file)
    filtered_structure.save(output_file, format='pqr')
    print(f"PQR file saved as '{output_file}'")
    
def split_pqr_file(pqr_file_path, hetatm_pqr_path, atom_pqr_path):
    """
    Splits a PQR file into separate HETATM and ATOM entries.

    Parameters:
    pqr_file_path (str): Path to the input PQR file.
    hetatm_pqr_path (str): Path to save the HETATM entries.
    atom_pqr_path (str): Path to save the ATOM entries.
    """
    # Read the contents of the PQR file
    with open(pqr_file_path, 'r') as pqr_file:
        pqr_contents = pqr_file.readlines()
    # Split the contents into HETATM and ATOM entries
    hetatm_lines = [line for line in pqr_contents if line.startswith("HETATM")]
    atom_lines = [line for line in pqr_contents if line.startswith("ATOM")]
    # Save HETATM entries to a separate PQR file
    with open(hetatm_pqr_path, 'w') as hetatm_file:
        hetatm_file.writelines(hetatm_lines)
    # Save ATOM entries to another PQR file
    with open(atom_pqr_path, 'w') as atom_file:
        atom_file.writelines(atom_lines)
    print(f"HETATM entries saved to: {hetatm_pqr_path}")
    print(f"ATOM entries saved to: {atom_pqr_path}")
    
def modify_ligand_pqr(file_path):
    """
    Reads a PQR file, modifies its content, and overwrites the original file.
    
    Args:
    file_path (str): Path to the PQR file.
    """
    try:
        # Read the original content
        with open(file_path, 'r') as file:
            pqr_lines = file.readlines()
        modified_lines = []
        for index, line in enumerate(pqr_lines):
            parts = line.split()
            # Update atom number and residue number
            parts[1] = str(index + 1)  # Atom number
            parts[5] = str(index + 1)  # Residue number
            # Remove chain ID and ensure correct formatting
            parts[4] = ''  # Chain ID is blank
            # Reformat the line with correct spacing
            modified_line = f"{parts[0]:<6}{parts[1]:>5}  {parts[2]:<4} {parts[3]:<4}{parts[4]:>1}{parts[5]:>4}  {parts[6]:>9}  {parts[7]:>6}  {parts[8]:>6} {parts[9]:>8} {parts[10]:>8}\n"
            modified_lines.append(modified_line)
        # Overwrite the original file with modified content
        with open(file_path, 'w') as file:
            file.writelines(modified_lines)
    except Exception as e:
        print(f"An error occurred: {e}")

def create_pdb_protein_ligand(input_pdb, output_pdb):
    """
    Filters out specified entries from a PDB file and saves the filtered content to a new file.

    Args:
        input_pdb (str): Path to the input PDB file.
        output_pdb (str): Path to the output PDB file.
    """
    entries = ["HOH", "WAT", "CRYST", "CONECT", "NA", "CL", "REMARK", "TER", "END"]
    # Read the input PDB file
    with open(input_pdb, 'r') as file:
        lines = file.readlines()
    # Filter lines to exclude certain entries
    filtered_lines = []
    for line in lines:
        if line.startswith(tuple(entries)):
            continue
        if line.startswith("ATOM") or line.startswith("HETATM"):
            residue_name = line[17:20].strip()
            if residue_name in entries:
                continue
        filtered_lines.append(line)
    # Write the filtered lines to the output PDB file
    with open(output_pdb, 'w') as file:
        file.writelines(filtered_lines)
    print(f"PDB file saved as '{output_pdb}'")

def extract_atoms_from_pdb(pdb_file):
    """
    Reads a PDB file and extracts atom types for protein and ligand atoms, counts each entry, and prints unique atom types.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        tuple: Two lists containing atom types for protein and ligand, respectively.
    """
    protein_atoms = []
    ligand_atoms = []
    # Read the PDB file
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                atom_type = line[12:16].strip()
                protein_atoms.append(atom_type[0] if atom_type[0].isalpha() else atom_type)
            elif line.startswith("HETATM"):
                atom_type = line[12:16].strip()
                ligand_atoms.append(atom_type[0] if atom_type[0].isalpha() else atom_type)
    # Create sets for unique atom types
    unique_protein_atoms = set(protein_atoms)
    unique_ligand_atoms = set(ligand_atoms)
    #print("Protein atom types:", protein_atoms)
    #print("Ligand atom types:", ligand_atoms)
    print("Total protein atoms:", len(protein_atoms))
    print("Total ligand atoms:", len(ligand_atoms))
    print("Unique protein atom types:", unique_protein_atoms)
    print("Unique ligand atom types:", unique_ligand_atoms)
    return protein_atoms, ligand_atoms

def map_atom_radii(atom_list, atom_radii):
    """
    Maps atom types to their corresponding radii.

    Args:
        atom_list (list): List of atom types.
        atom_radii (dict): Dictionary of atom types to their corresponding radii.

    Returns:
        list: List of radii corresponding to the atom types in the input list.
    """
    radii_list = []
    for atom in atom_list:
        if atom in atom_radii:
            radii_list.append(atom_radii[atom])
        else:
            radii_list.append(None)  
    return radii_list

def update_pqr_radii_inplace(pqr_file_path, new_radii_list):
    """
    Updates the radii in a PQR file with new radii values, overwriting the original file while preserving original formatting.

    Args:
        pqr_file_path (str): Path to the PQR file to be updated.
        new_radii_list (list): List of new radii values.
    """
    updated_lines = []
    with open(pqr_file_path, 'r') as file:
        lines = file.readlines()
    for line, new_radii in zip(lines, new_radii_list):
        # Locate the last whitespace before the radius and slice the string there
        last_space_index = line.rfind(' ') + 1
        # Preserve the formatting by reconstructing the line with the new radius
        updated_line = line[:last_space_index] + f"{new_radii:.4f}\n"
        updated_lines.append(updated_line)
    # Overwrite the original file with the updated radii
    with open(pqr_file_path, 'w') as file:
        file.writelines(updated_lines)
    print(f"Updated PQR file saved over the original file at {pqr_file_path}")

# Function: Load system from an XML file, coordinates from a PDB file, filter out unwanted residues, and save the filtered structure as a PQR file.
# 1. Loads a molecular system from an XML file, which typically contains force field information and system configuration.
# 2. Loads coordinate information from a PDB file, which provides the positions of atoms.
# 3. Filters out unwanted residues such as water (HOH, WAT) and ions (NA, CL) to simplify the system for certain types of analyses or simulations.
# 4. Saves the filtered system in a PQR format, which includes atomic charges and radii along with the coordinates, commonly used in electrostatics calculations.
# 5. The output is stored in 'complex.pqr', effectively preparing it for further processing or simulation tasks.
get_pqr(xml_file='complex_serialized.xml', pdb_file='complex_SEEKR.pdb', output_file='complex.pqr')

# Function: Splits a PQR file into separate files for HETATM and ATOM entries.
# 1. Reads the 'complex.pqr' file which contains both protein and ligand information.
# 2. Extracts lines starting with "HETATM", typically representing ligand or non-standard residues, and saves them to 'ligand.pqr'.
# 3. Extracts lines starting with "ATOM", typically representing standard protein residues, and saves them to 'receptor.pqr'.
# 4. This separation facilitates focused analyses or simulations on either the protein or the ligand.
split_pqr_file("complex.pqr", "ligand.pqr", "receptor.pqr")

# Function: Modifies a PQR file to renumber atom and residue indices and removes chain IDs, then saves the changes.
# 1. Opens 'ligand.pqr', which contains ligand information possibly from a larger complex.
# 2. Reads and modifies each line to update atom and residue numbers sequentially, ensuring each entry is unique and sequential from 1 onwards.
# 3. Removes chain IDs to simplify the data, possibly for use in simulation setups where chain distinction is unnecessary.
# 4. Rewrites the file with the modified content, ensuring the data format is consistent and aligned with the requirements of most molecular simulation software.
# 5. This update is critical for ensuring that simulations or analyses using this file do not encounter errors due to improper indexing or format discrepancies.
modify_ligand_pqr(file_path="ligand.pqr")

# Remove unwanted residues and entries from a PDB file and save the cleaned version.
# 1. Reads an input PDB file ('complex_minimized.pdb') containing molecular structures, possibly including unwanted residues like water and ions.
# 2. Filters out specified entries such as water molecules, ions, and annotation lines that are not needed for subsequent analyses.
# 3. Saves the cleaned and streamlined molecular structure to a new PDB file ('receptor_ligand.pdb'), facilitating more efficient processing in later steps.
create_pdb_protein_ligand(input_pdb='complex_SEEKR.pdb', output_pdb='receptor_ligand.pdb')

# Extracts atom types from the cleaned PDB file, separating them into protein and ligand categories.
# 1. Opens the previously cleaned PDB file to read the structural data.
# 2. Iterates through each line, classifying atom entries into protein (ATOM) and ligand (HETATM) categories based on their line prefixes.
# 3. Extracts and lists the types of atoms found in each category, providing essential data for further molecular characterization and analysis.
protein_atoms, ligand_atoms = extract_atoms_from_pdb("receptor_ligand.pdb")

# Dictionary of atom types mapped to their corresponding radii, typically used in molecular modeling.
# 1. Establishes a mapping of common atom types to their Van der Waals radii, which are crucial for simulations involving molecular interactions and dynamics.
# 2. This dictionary will be used to apply these radii values to the respective atom types in the structure files, aligning with physical and chemical properties.
atom_radii = {'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80, 'F': 1.47, 'Cl': 1.75, 
              'Br': 1.85, 'I': 1.98, 'Fe': 1.80, 'Cu': 1.40, 'Zn': 1.39, 'Ca': 1.00, 'Mg': 1.18}

# Map the extracted atom types to their corresponding radii for protein and ligand atoms.
# 1. Converts lists of atom types extracted from the PDB file into lists of corresponding radii using the predefined dictionary.
# 2. Each atom type in the protein and ligand lists is associated with its radius, essential for accurate representation in simulations and visualizations.
protein_atom_radii = map_atom_radii(atom_list=protein_atoms, atom_radii=atom_radii)
ligand_atom_radii = map_atom_radii(atom_list=ligand_atoms, atom_radii=atom_radii)

# Print the mapped radii lists for proteins and ligands for verification.:
# 1. Outputs the computed radii lists for both proteins and ligands, allowing for a quick verification and validation of the mapping process.
print("Protein Atom Radii:", protein_atom_radii)  
print("Ligand Atom Radii:", ligand_atom_radii)  

# Update the receptor and ligand PQR file in place with the new radii, maintaining the original file format.
# 1. Reads the PQR files for receptor and ligand, which initially contain standard radii values.
# 2. Updates each line with new radii values from the mapped lists, ensuring that the physical properties represented are accurate and tailored to the specific atoms.
# 3. Saves the modifications directly over the original files, preserving the format but updating the content to reflect the new radii, essential for accurate computational assessments.
update_pqr_radii_inplace(pqr_file_path='receptor.pqr', new_radii_list=protein_atom_radii)
update_pqr_radii_inplace(pqr_file_path='ligand.pqr', new_radii_list=ligand_atom_radii)

# Remove any intermediate files
