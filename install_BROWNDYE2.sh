#!/bin/bash

# Step 1: Navigate to the SageMaker directory
cd /home/ec2-user/SageMaker

# Step 2: Install the EPEL repository
sudo amazon-linux-extras install epel -y

# Step 3: Install necessary development tools and dependencies
sudo yum install -y devtoolset-9 ocaml expat-devel lapack-devel apbs

# Step 4: If the above command fails, try with --skip-broken to avoid broken packages
if [ $? -ne 0 ]; then
    sudo yum install -y devtoolset-9 ocaml expat-devel lapack-devel apbs --skip-broken
fi

# Step 5: Download the Browndye2 package
rm -rf browndye2.tar.gz
rm -rf browndye2
wget https://browndye.ucsd.edu/downloads/browndye2.tar.gz

# Step 6: Extract the package
tar -xvzf browndye2.tar.gz

# Step 7: Remove the tar.gz file to save space
rm -rf browndye2.tar.gz

# Step 8: Navigate into the Browndye2 directory
cd browndye2

# Step 9: Compile the Browndye2 software using 4 cores
make -j 4 all

# Step 10: Return to the SageMaker directory
cd /home/ec2-user/SageMaker

