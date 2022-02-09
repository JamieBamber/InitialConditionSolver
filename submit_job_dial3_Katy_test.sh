#!/bin/bash
#
# This script should make a new directory in GRChombo_data for the data, and copy the slurm_submit file, params.txt file and output files to
# that directory, submit the job, then change back to the working directory

# this copy is for COSMA6

work_dir=/home/dc-bamb1/InitialConditionSolver
cd $work_dir
data_directory=/scratch/dp092/dc-bamb1/GRChombo_data/InitialConditionSolver

new_dir=Katys_test
new_dir_path=${data_directory}/${new_dir}
#
mkdir -p ${new_dir_path}

cd ${work_dir}
cp pbs_submit_dial3 ${new_dir_path}/pbs_submit

#params_file=params_ratio${ratio}.txt
params_file=params.txt

cp ${params_file} ${new_dir_path}/params.txt

cd ${new_dir_path}
mkdir -p outputs
cd outputs
sbatch ../pbs_submit
#
cd ${work_dir}
