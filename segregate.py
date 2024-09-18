from Bio import PDB

input_pdb = "test.pdb"
protein_pdb = "protein.pdb"
ligand_res_name = "INH" 
ligand_pdb = "ligand.pdb"

# Remove ions and water, keeping only the protein and ligand
parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure("protein_ligand", input_pdb)

class NonWaterNonIonSelect(PDB.Select):
    """ Selects only protein atoms and excludes water and ions. """
    def accept_residue(self, residue):
        # Exclude water and common ions (identified by hetero-flag or water ID)
        return residue.id[0] == " " and residue.resname != ligand_res_name  # Accept protein residues only

class LigandSelect(PDB.Select):
    """ Selects only the ligand based on its residue name. """
    def accept_residue(self, residue):
        return residue.resname == ligand_res_name  # Select only the ligand residue

# Save the protein to a new PDB file
io = PDB.PDBIO()
io.set_structure(structure)
io.save(protein_pdb, NonWaterNonIonSelect())

# Extract ligand and ensure correct valencies using OpenEye

# Save ligand to a temporary PDB file
io.save(ligand_pdb, LigandSelect())
