import xml.etree.ElementTree as ET
import os

# Determine the current working directory
current_directory = os.getcwd()

# Define the path to the SEEKR_SIMULATION folder
seekr_simulation_dir = os.path.join(current_directory, "SEEKR_SIMULATION")

# Ensure the SEEKR_SIMULATION directory exists
if not os.path.exists(seekr_simulation_dir):
    print(f"The directory {seekr_simulation_dir} does not exist.")
    exit(1)

# Define the path to the model.xml file
model_xml_path = os.path.join(seekr_simulation_dir, "model.xml")

# Load the original XML file
tree = ET.parse(model_xml_path)
root = tree.getroot()

# Find the calculation settings element
calculation_settings = root.find(".//calculation_settings")

# Modify the specified elements
if calculation_settings is not None:
    calculation_settings.find(".//num_production_steps").text = "200000"
    calculation_settings.find(".//energy_reporter_interval").text = "10000"
    calculation_settings.find(".//trajectory_reporter_interval").text = "10000"
    calculation_settings.find(".//restart_checkpoint_interval").text = "10000"

# Write the modified XML to a new file in the SEEKR_SIMULATION directory
sample_model_path = os.path.join(seekr_simulation_dir, "sample_model.xml")
tree.write(sample_model_path)

print(f"The file 'sample_model.xml' has been created in {seekr_simulation_dir} with the updated settings.")

