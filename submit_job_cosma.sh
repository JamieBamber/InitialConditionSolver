#!/bin/bash
#
# This script should make a new directory in GRChombo_data for the data, and copy the slurm_submit file, params.txt file and output files to
# that directory, submit the job, then change back to the working directory

# this copy is for COSMA6

work_dir=/cosma/home/dp174/dc-bamb1/InitialConditionSolver_complex
cd $work_dir
data_directory=/cosma6/data/dp174/dc-bamb1/GRChombo_data/InitialConditionSolver

run0001=(0.48847892320123 12.21358 0.34 0.1356676906 0 0 0)

#params_file=params_ratio${ratio}.txt
params_file=params.txt

G=1 # 10^{-10}

run_list=(
	run0001
)

r_inner=2
L=512
N1=64

for run in "${run_list[@]}"
do
  	cd $work_dir
        # extract parameters    
        val="$run[0]"; M="${!val}"
        val="$run[1]"; d="${!val}"
        val="$run[2]"; mu="${!val}"
        val="$run[3]"; dt_mult="${!val}"
        val="$run[4]"; l="${!val}"
        val="$run[5]"; m="${!val}"
        val="$run[6]"; Al="${!val}"

        #omega_BH=0
        omega_BH=$(awk "BEGIN {printf \"%.7f\n\", sqrt(2*${M}/(${d}*${d}*${d}))}")
        #omega_BH=0.25
        echo "omega_BH = ${omega_BH}"

	new_dir=gaussian_kappa0.0125_mu${mu}_l0_m0_G${G}_max_level9
	#Homogeneous_${run}_mu${mu}_l0_m0_G${G}_max_level9
	#Newtonian_${run}_G${G}_max_level9_wslope0.25_wradius50_n${num}
	new_dir_path=${data_directory}/${new_dir}
	#
	mkdir -p ${new_dir_path}
	echo "made ${new_dir_path}"

	cd ${work_dir}
	cp slurm_submit_cosma ${new_dir_path}/slurm_submit
	
	cp ${params_file} ${new_dir_path}/params.txt

	#sed -i "s|SUBDIR|${subdir}|" ${new_dir_path}/params.txt
	sed -i "s|GVAL|${G}|" ${new_dir_path}/params.txt
	sed -i "s|OMEGAB|${omega_BH}|" ${new_dir_path}/params.txt
	sed -i "s|OUTNAME|${new_dir}|" ${new_dir_path}/params.txt
	sed -i "s|MUVAL|${mu}|" ${new_dir_path}/params.txt
	sed -i "s|NPLOTFL|${num}|" ${new_dir_path}/params.txt

	cd ${new_dir_path}
	mkdir -p outputs
        cd outputs

	sbatch ../slurm_submit
	#
	cd ${work_dir}
done
