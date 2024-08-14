#!/bin/bash

# Step 1: Change to the SageMaker directory
cd /home/ec2-user/SageMaker

# Step 2: Remove existing environments and directories
rm -rf /home/ec2-user/SageMaker/envs/SEEKR2	
rm -rf /home/ec2-user/SageMaker/seekr2
rm -rf /home/ec2-user/SageMaker/seekrtools

# Step 3: Source the bashrc file
source /home/ec2-user/SageMaker/.bashrc

# Step 4: Remove the SEEKR2 environment if it still exists
conda remove --prefix /home/ec2-user/SageMaker/envs/SEEKR2 --all --yes

# Step 5: Create a new SEEKR2 environment
conda create --prefix /home/ec2-user/SageMaker/envs/SEEKR2 python=3.9 --yes

# Step 6: Activate the new environment
source activate /home/ec2-user/SageMaker/envs/SEEKR2

# Step 7: Install required Python packages with pip
/home/ec2-user/SageMaker/envs/SEEKR2/bin/pip install cython
/home/ec2-user/SageMaker/envs/SEEKR2/bin/pip install mpi4py
/home/ec2-user/SageMaker/envs/SEEKR2/bin/pip install swig

# Step 8: Install required Conda packages
conda install conda-forge::git --yes
conda install conda-forge::netcdf4 --yes
conda install conda-forge::scipy --yes
conda install conda-forge::numpy --yes
conda install conda-forge::doxygen --yes
conda install conda-forge::mdtraj --yes
conda install conda-forge::ambertools --yes

# Step 9: Install the SEEKR2 OpenMM plugin
conda install -c conda-forge seekr2_openmm_plugin --yes

# Step 10: Verify installation
python -c "import seekr2plugin"

# Step 11: Clone the SEEKR2 repository and install
cd /home/ec2-user/SageMaker
git clone https://github.com/seekrcentral/seekr2.git
cd /home/ec2-user/SageMaker/seekr2
/home/ec2-user/SageMaker/envs/SEEKR2/bin/pip install .

# Step 12: Run tests for SEEKR2
pytest

# Step 13: Clone the SEEKRTools repository and install
cd /home/ec2-user/SageMaker
git clone https://github.com/seekrcentral/seekrtools.git
cd /home/ec2-user/SageMaker/seekrtools
/home/ec2-user/SageMaker/envs/SEEKR2/bin/pip install .

# Step 14: Run tests for SEEKRTools
pytest

# Step 15: Change to the SageMaker  directory
cd /home/ec2-user/SageMaker
