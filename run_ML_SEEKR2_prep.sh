#!/bin/bash

# Export the OpenEye license
export OE_LICENSE="/home/aojha/licenses/oe_license.txt"
echo "OpenEye license has been exported."

# Set explicit paths for SEEKR2 and SEEKRTOOLS directories
SEEKR2_DIR="/home/aojha/seekr2"
SEEKRTOOLS_DIR="/home/aojha/seekrtools"

# Exit immediately if a command exits with a non-zero status
set -e

# Initialize Conda (adjust the path to the conda.sh script if needed)
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Activate the ESPALOMA environment
echo "Activating ESPALOMA environment..."
conda activate ESPALOMA

# Run parameterization script
echo "Running parameterize.py..."
python parameterize.py

# Run write input XML script
echo "Running write_input_xml.py..."
python write_input_xml.py

# Run generate Browndye files script
echo "Running generate_browndye_files.py..."
python generate_browndye_files.py

# Activate the SEEKR2 environment
echo "Activating SEEKR2 environment..."
conda activate SEEKR2

# Check if the SEEKR2 and SEEKRTOOLS directories exist
if [[ ! -d "$SEEKR2_DIR" ]]; then
  echo "Error: SEEKR2 directory not found at $SEEKR2_DIR."
  exit 1
fi

if [[ ! -d "$SEEKRTOOLS_DIR" ]]; then
  echo "Error: SEEKRTOOLS directory not found at $SEEKRTOOLS_DIR."
  exit 1
fi

echo "Using SEEKR2 directory at: $SEEKR2_DIR"
echo "Using SEEKRTOOLS directory at: $SEEKRTOOLS_DIR"

# Run SEEKR2 prepare step using the explicitly set path
echo "Running SEEKR2 prepare.py..."
python "$SEEKR2_DIR/seekr2/prepare.py" input.xml

# Run HIDR simulation with SEEKR2 using the explicitly set path
echo "Running HIDR simulation..."
python "$SEEKRTOOLS_DIR/seekrtools/hidr/hidr.py" any SEEKR_SIMULATION/model.xml -M metaD -p complex_SEEKR.pdb -b 5

# Create the sample_model.xml file by running create_sample_model.py
echo "Creating sample_model.xml..."
python create_sample_model.py

# Run SEEKR2 simulation with sample_model.xml
echo "Running SEEKR2 simulation with sample_model.xml..."
python "$SEEKR2_DIR/seekr2/run.py" any SEEKR_SIMULATION/sample_model.xml

# Delete the sample_model.xml file after the simulation
echo "Deleting sample_model.xml..."
rm SEEKR_SIMULATION/sample_model.xml

echo "All commands executed successfully."

