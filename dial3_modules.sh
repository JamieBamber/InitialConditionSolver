#!/bin/bash

module purge

# Load your desired modules here
module load arm/forge/21.0.2
module load cray-pmi cray-pmi-lib
#module swap cray-mpich/8.1.4 cray-mpich-ucx/8.1.4
#module swap craype-network-ofi craype-network-ucx
#module load craype-network-infiniband
#module load craype-network-ofi
module load craype
module load craype-network-ucx
module load cce/11.0.4
module load cray-mpich-ucx/8.1.4

module load cray-libsci/21.04.1.1
#module load cray-hdf5-parallel/1.12.0.2

#module load gcc/11.1.0/xu5zmz
module load cray-hdf5/1.12.0.2
module load intel-oneapi-compilers/2021.2.0/3q6eev
module load intel-parallel-studio/cluster.2019.5/aveeft
module load cray-cti/2.12.2


