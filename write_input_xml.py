import numpy as np
import os

def process_pdb(input_pdb, output_pdb):
    # Keywords to filter out
    keywords = ["REMARK", "CRYST", "TER", "HOH", "WAT", "CONECT", "END", "Na", "Cl" ]
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if not any(keyword in line for keyword in keywords):
                outfile.write(line)

def get_full_path(filename):
    # Get the current working directory
    current_directory = os.getcwd()
    # Join the current directory with the filename
    full_path = os.path.join(current_directory, filename)
    return full_path

def parse_pdb(filename):
    """Parse the PDB file to separate protein and ligand atoms."""
    with open(filename, 'r') as file:
        lines = file.readlines()
    protein_atoms = []
    ligand_atoms = []
    for line in lines:
        if line.startswith("HETATM") and " CA " in line:  
            atom_index = int(line[6:11].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            protein_atoms.append((atom_index, np.array([x, y, z])))
        elif line.startswith("HETATM") and "XX1" in line:
            atom_index = int(line[6:11].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            ligand_atoms.append((atom_index, np.array([x, y, z])))
    return protein_atoms, ligand_atoms

def extract_alpha_carbon_and_ligand_indices(protein_atoms, ligand_atoms, threshold):
    """Extract indices of alpha carbons close to the ligand and indices of all ligand atoms."""
    close_c_alpha_atom_indices = []
    ligand_atom_indices = [atom_index for atom_index, _ in ligand_atoms]
    for atom_index, ca_pos in protein_atoms:
        for _, lig_pos in ligand_atoms:
            distance = np.linalg.norm(ca_pos - lig_pos)
            if distance <= threshold:
                close_c_alpha_atom_indices.append(atom_index)
                break  
    return close_c_alpha_atom_indices, ligand_atom_indices

def adjust_alpha_carbon_indices(alpha_carbon_indices):
    """Adjust alpha carbon indices by subtracting one and calculate the length."""
    adjusted_indices = [index - 1 for index in alpha_carbon_indices]
    length_of_indices = len(adjusted_indices)
    return adjusted_indices, length_of_indices

def adjust_ligand_indices(ligand_indices):
    """Adjust ligand indices by subtracting one and calculate the length."""
    adjusted_indices = [index - 2 for index in ligand_indices]
    length_of_indices = len(adjusted_indices)
    return adjusted_indices, length_of_indices

def create_sequential_ligand_list(length):
    """Create a sequential list of integers from 0 to length-1."""
    return list(range(length))

def write_input_xml(filename,
                    calculation_type,
                    md_output_interval,
                    md_steps_per_anchor,
                    temperature,
                    pressure,
                    ensemble,
                    root_directory,
                    md_program,
                    constraints,
                    rigidWater,
                    hydrogenMass,
                    timestep,
                    nonbonded_cutoff,
                    receptor_indices,
                    ligand_indices_openMM,
                    radii,
                    system_filename,
                    receptor_pqr_filename,
                    ligand_pqr_filename,
                    ligand_indices_BD,
                    num_b_surface_trajectories,
                    n_threads,
                    comment_browndye_settings=False):  # Added parameter to control commenting

    # Convert the lists to strings for XML formatting
    receptor_indices_str = ", ".join(map(str, receptor_indices))
    ligand_indices_openMM_str = ", ".join(map(str, ligand_indices_openMM))
    ligand_indices_BD_str = ", ".join(map(str, ligand_indices_BD))

    # Initialize the input_anchors content
    input_anchors_content = ""
    # Loop through the radii to generate the XML content
    for radius in radii:
        bulk_anchor = "False"
        # Check if it is the last radius to set bulk_anchor to True
        if radius == radii[-1]:
            bulk_anchor = "True"
            system_filename = ""  # Assuming the last one does not have a system_filename
        # Append the input_anchor content for each radius
        input_anchors_content += f"""
                <input_anchor class="Spherical_cv_anchor">
                    <radius>{radius}</radius>
                    <lower_milestone_radius/>
                    <upper_milestone_radius/>
                    <starting_forcefield_params class="Forcefield_params">
                        <system_filename>{system_filename}</system_filename>
                        <box_vectors/>
                        <pdb_coordinates_filename></pdb_coordinates_filename>
                    </starting_forcefield_params>
                    <bound_state>False</bound_state>
                    <bulk_anchor>{bulk_anchor}</bulk_anchor>
                </input_anchor>"""

    # Handle the hydrogenMass tag based on the value
    if hydrogenMass == 1:
        hydrogenMass_tag = "<hydrogenMass/>"
    else:
        hydrogenMass_tag = f"<hydrogenMass>{hydrogenMass}</hydrogenMass>"

    # Construct ions content
    ion_details = [
        {"radius": 1.2, "charge": -1.0, "conc": 0.01},
        {"radius": 0.9, "charge": 1.0, "conc": 0.01}
    ]
    ions_content = ""
    for ion in ion_details:
        ions_content += f"""
            <ion class="Ion">
                <radius>{ion['radius']}</radius>
                <charge>{ion['charge']}</charge>
                <conc>{ion['conc']}</conc>
            </ion>"""

    # Define the browndye_settings_input content
    browndye_settings_input_content = f"""
    <browndye_settings_input class="Browndye_settings_input">
        <binary_directory></binary_directory>
        <receptor_pqr_filename>{receptor_pqr_filename}</receptor_pqr_filename>
        <ligand_pqr_filename>{ligand_pqr_filename}</ligand_pqr_filename>
        <apbs_grid_spacing>0.5</apbs_grid_spacing>
        <receptor_indices>[{receptor_indices_str}]</receptor_indices>
        <ligand_indices>[{ligand_indices_BD_str}]</ligand_indices>
        <ions>{ions_content}
        </ions>
        <num_b_surface_trajectories>{num_b_surface_trajectories}</num_b_surface_trajectories>
        <n_threads>{n_threads}</n_threads>
    </browndye_settings_input>"""

    # Conditionally comment out the browndye_settings_input section
    if comment_browndye_settings:
        browndye_settings_input_content = f"<!--{browndye_settings_input_content}-->"

    # Combine all parts into the final XML content
    xml_content = f"""<?xml version="1.0" ?>
<model_input class='Model_input'>
    <calculation_type>{calculation_type}</calculation_type>
    <calculation_settings class="MMVT_input_settings">
        <md_output_interval>{md_output_interval}</md_output_interval>
        <md_steps_per_anchor>{md_steps_per_anchor}</md_steps_per_anchor>
    </calculation_settings>
    <temperature>{temperature}</temperature>
    <pressure>{pressure}</pressure>
    <ensemble>{ensemble}</ensemble>
    <root_directory>{root_directory}</root_directory>
    <md_program>{md_program}</md_program> 
    <constraints>{constraints}</constraints>
    <rigidWater>{rigidWater}</rigidWater>
    {hydrogenMass_tag}
    <timestep>{timestep}</timestep>
    <nonbonded_cutoff>{nonbonded_cutoff}</nonbonded_cutoff>
    <cv_inputs>
        <cv_input class="Spherical_cv_input">
            <group1>[{receptor_indices_str}]</group1>
            <group2>[{ligand_indices_openMM_str}]</group2>
            <input_anchors>{input_anchors_content}
            </input_anchors>
        </cv_input>
    </cv_inputs>
    {browndye_settings_input_content}
</model_input>"""

    # Write the final content to the input.xml file
    with open(filename, "w") as f:
        f.write(xml_content)

process_pdb(input_pdb="complex_minimized.pdb", output_pdb="receptor_ligand.pdb")
# Parsing the PDB file to get protein and ligand atoms
protein_atoms, ligand_atoms = parse_pdb(filename='receptor_ligand.pdb')
# Extracting receptor and ligand indices
receptor_alpha_indices, ligand_indices = extract_alpha_carbon_and_ligand_indices(protein_atoms=protein_atoms,ligand_atoms=ligand_atoms, threshold=6.00)
adjusted_receptor_alpha_indices, receptor_alpha_carbon_length = adjust_alpha_carbon_indices(receptor_alpha_indices)
print("Receptor Alpha Carbon Indices:", receptor_alpha_indices)
print("Receptor Alpha Carbon Indices for model.xml SEEKR input file:", adjusted_receptor_alpha_indices)
print("Number of Receptor Alpha-Carbon Atoms:", receptor_alpha_carbon_length)
adjusted_ligand_indices, ligand_length = adjust_ligand_indices(ligand_indices)
browndye_ligand_list = create_sequential_ligand_list(ligand_length)
print("Ligand Atom Indices:", ligand_indices)
print("Ligand Atom Indices for model.xml SEEKR input file:", adjusted_ligand_indices)
print("Number of Ligand Atoms:", ligand_length)
print("Browndye2 Ligand Indices:", browndye_ligand_list)

write_input_xml(
    filename="input.xml",
    calculation_type="mmvt",
    md_output_interval=100000,
    md_steps_per_anchor=100000000,
    temperature=300.0,
    pressure=1.0,
    ensemble="nvt",
    root_directory=get_full_path(filename="SEEKR_SIMULATION"),
    md_program="openmm",
    constraints="HBonds",
    rigidWater="True",
    hydrogenMass=1,
    timestep=0.002,
    nonbonded_cutoff=1.0,
    receptor_indices=adjusted_receptor_alpha_indices,
    ligand_indices_openMM=adjusted_ligand_indices,
    radii=[0.125, 0.2, 0.275, 0.35, 0.425, 0.5, 0.575, 0.65, 0.725, 0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0],
    system_filename=get_full_path(filename="complex_serialized.xml"),
    receptor_pqr_filename=get_full_path(filename="receptor.pqr"),
    ligand_pqr_filename=get_full_path(filename="ligand.pqr"),
    ligand_indices_BD=browndye_ligand_list,
    num_b_surface_trajectories=100000,
    n_threads=12,
    comment_browndye_settings=True 
)


