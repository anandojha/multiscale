#!/bin/bash

# Remove existing SEEKR2 conda environment
conda env remove --name SEEKR2 --yes

# Navigate to the home directory
cd ~

# Remove existing seekr2 and seekrtools directories
rm -rf seekr2 seekrtools

# Create a new conda environment with Python 3.9
conda create --name SEEKR2 python=3.9 --yes

# Activate the SEEKR2 environment
source activate SEEKR2

# Upgrade Cython using pip
pip install --upgrade cython

# Install packages with conda
conda install git numpy scipy netcdf4 mpi4py swig --yes
conda install -c conda-forge ambertools mdtraj doxygen --yes

# Install cmake curses gui
sudo apt-get install cmake-curses-gui

# Install seekr2_openmm_plugin from conda-forge
conda install -c conda-forge seekr2_openmm_plugin --yes

# Clone the seekr2 repository
git clone https://github.com/seekrcentral/seekr2.git

# Install seekr2
cd seekr2
python -m pip install .
# Run tests
pytest
cd ..

# Clone the seekrtools repository
git clone https://github.com/seekrcentral/seekrtools.git
cd seekrtools
# Install seekrtools
python -m pip install .
# Run tests
pytest
cd ..

