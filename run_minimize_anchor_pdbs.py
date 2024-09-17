import subprocess
import os

def remove_prod_files(seekr_simulation_dir):
    """
    Remove all files in anchor_*/prod/* directories.
    """
    try:
        # Use subprocess to run the shell command to remove files
        command = "rm anchor_*/prod/*"
        subprocess.run(command, shell=True, check=True, cwd=seekr_simulation_dir)
        print("Successfully removed all files in anchor_*/prod/* directories.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while removing files: {e}")

def run_minimize_anchor_pdbs(base_directory, script_name):
    # Store the current directory to return to it later
    current_directory = os.getcwd()
    # Step 1: Identify the SEEKR_SIMULATION directory
    seekr_simulation_dir = os.path.join(base_directory, "SEEKR_SIMULATION")
    if not os.path.exists(seekr_simulation_dir):
        raise FileNotFoundError(f"SEEKR_SIMULATION directory not found in {base_directory}")
    # Step 2: Identify all 'anchor' directories
    anchor_dirs = [d for d in os.listdir(seekr_simulation_dir) 
                   if os.path.isdir(os.path.join(seekr_simulation_dir, d)) and d.startswith('anchor_') and d[7:].isdigit()]
    # Step 3: Sort anchor directories based on the numeric part after 'anchor_'
    anchor_dirs.sort(key=lambda x: int(x[7:]))  # Sorting by the number after 'anchor_'
    # Step 4: Iterate over each 'anchor' directory
    for anchor_dir in anchor_dirs:
        anchor_path = os.path.join(seekr_simulation_dir, anchor_dir)
        building_path = os.path.join(anchor_path, "building")
        # Check if 'building' directory exists
        if not os.path.exists(building_path):
            print(f"Building directory not found in {anchor_dir}, skipping...")
            continue
        print(f"Processing {building_path}...")
        # Change the current working directory to the 'building' directory
        os.chdir(building_path)
        # Step 5: Run the minimize_anchor_pdbs.py script
        try:
            # Use subprocess to run the script located in the original directory
            result = subprocess.run(['python', os.path.join(current_directory, script_name)], capture_output=True, text=True)
            print(result.stdout)
            if result.stderr:
                print(f"Errors in {building_path}:\n{result.stderr}")
        except Exception as e:
            print(f"An error occurred while processing {building_path}: {e}")
        finally:
            # Change back to the original working directory after processing
            os.chdir(current_directory)

if __name__ == "__main__":
    # Get the current working directory
    base_directory = os.getcwd()
    # Specify the name of the Python script to run
    script_name = "minimize_anchor_pdbs.py"
    # Run the script in all identified anchor folders
    run_minimize_anchor_pdbs(base_directory, script_name)
