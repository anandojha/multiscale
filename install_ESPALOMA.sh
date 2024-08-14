#!/bin/bash

# Step 1: Change to the envs directory
cd /home/ec2-user/SageMaker/envs

# Step 2: Remove existing ESPALOMA environment
rm -rf /home/ec2-user/SageMaker/envs/ESPALOMA

# Step 3: Source the bashrc file
source /home/ec2-user/SageMaker/.bashrc

# Step 4: Remove the ESPALOMA environment if it still exists
conda remove --prefix /home/ec2-user/SageMaker/envs/ESPALOMA --all --yes

# Step 5: Create the espaloma_env.yml file
cat <<EOT > espaloma_env.yml
name: ESPALOMA
channels:
  - conda-forge
  - openeye
dependencies:
  - python=3.10
  - openff-toolkit
  - openmmforcefields
  - openeye-toolkits
  - espaloma=0.3.2
EOT

# Step 6: Create the Conda environment with Python 3.10
conda create --prefix /home/ec2-user/SageMaker/envs/ESPALOMA python=3.10 -y

# Step 7: Source the bashrc to ensure Conda is available
source /home/ec2-user/SageMaker/.bashrc

# Step 8: Activate the new environment
conda activate /home/ec2-user/SageMaker/envs/ESPALOMA

# Step 9: Install Mamba in the new environment
conda install -n base  mamba -y

# Step 10: Use Mamba to update the environment with the packages specified in the YAML file
mamba env update -f espaloma_env.yml -p /home/ec2-user/SageMaker/envs/ESPALOMA

# Step 11: Set the OpenEye license environment variable
export OE_LICENSE="/home/ec2-user/SageMaker/oe_license.txt"
echo "Environment setup complete and OE_LICENSE set."

# Step 12: Source the bashrc to ensure conda is available
source /home/ec2-user/SageMaker/.bashrc

# Step 13: Activate the new environment
conda activate /home/ec2-user/SageMaker/envs/ESPALOMA

# Step 14: Remove the espaloma_env.yml file
rm -rf espaloma_env.yml

# Step 15: Check if the required packages can be imported successfully
python - <<EOF
try:
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
    print("All imports successful and ready to run ESPALOMA")
except ImportError as e:
    print(f"ImportError: {e}")
    exit(1)
EOF
