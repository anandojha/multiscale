from Bio.PDB import PDBParser
from rdkit import Chem
import numpy as np

def calculate_com(coords, masses):
    total_mass = np.sum(masses)
    com = np.sum(coords.T * masses, axis=1) / total_mass
    return com

def get_sdf_com(sdf_file):
    suppl = Chem.SDMolSupplier(sdf_file)
    mol = suppl[0]  # Assuming the first molecule

    coords = mol.GetConformer().GetPositions()
    atom_weights = [atom.GetMass() for atom in mol.GetAtoms()]

    return calculate_com(np.array(coords), np.array(atom_weights))

def get_pdb_com(pdb_file, alpha_carbon_atom_indices):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)

    alpha_carbon_coords = []
    alpha_carbon_masses = []

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_id() == 'CA' and atom.get_serial_number() in alpha_carbon_atom_indices:
                        alpha_carbon_coords.append(atom.get_coord())
                        alpha_carbon_masses.append(12.01)  # Approximate mass of carbon

    return calculate_com(np.array(alpha_carbon_coords), np.array(alpha_carbon_masses))

def com_com_distance(sdf_com, pdb_com):
    return np.linalg.norm(sdf_com - pdb_com)

# Path to your SDF and PDB files
sdf_file_path = "ligand.sdf"
pdb_file_path = "protein.pdb"

# Receptor Alpha Carbon Indices
alpha_carbon_atom_indices = [530, 579, 593, 615, 627, 637, 668, 1253, 1265, 1279, 1286, 1305, 1312, 1329, 
                             1375, 1387, 1439, 1453, 1472, 1479, 1862, 1878, 1885, 1905, 1937, 2050, 2066, 
                             2080, 2096, 2115, 2129, 2151, 2259, 2269, 2563, 2570, 2584, 2606]

# Get COM for the SDF file atoms
sdf_com = get_sdf_com(sdf_file_path)

# Get COM for selected alpha carbon atoms in the PDB file
pdb_com = get_pdb_com(pdb_file_path, alpha_carbon_atom_indices)

# Calculate the COM-COM distance
com_distance = com_com_distance(sdf_com, pdb_com)

print(f"The COM-COM distance is: {com_distance} Å")

