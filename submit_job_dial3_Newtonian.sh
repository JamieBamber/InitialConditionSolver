#!/bin/bash
#
# This script should make a new directory in GRChombo_data for the data, and copy the slurm_submit file, params.txt file and output files to
# that directory, submit the job, then change back to the working directory

# this copy is for COSMA6

work_dir=/home/dc-bamb1/InitialConditionSolver_complex
cd $work_dir
data_directory=/scratch/dp092/dc-bamb1/GRChombo_data/InitialConditionSolver

run0031=(0.48847892320123 12.21358 0.34 0.1356676906 0 0 0)
run0032=(0.48847892320123 12.21358 0.34 0.1356676906 1 1 0)
run0033=(0.48847892320123 12.21358 0.34 0.1356676906 1 -1 0)

run0034=(0.48847892320123 12.21358 0.5 0.1356676906 0 0 0)

run0035=(0.48847892320123 12.21358 0.01 0.2713353812 0 0 0)

run0037=(0.48847892320123 12.21358 1.0 0.0678338453 0 0 0)
run0038=(0.48847892320123 12.21358 1.0 0.0678338453 1 1 0)
run0039=(0.48847892320123 12.21358 1.0 0.0678338453 1 -1 0)
run0040=(0.48847892320123 0 1.0 0.0678338453 0 0 0)
run0041=(0.48847892320123 0 1.0 0.0678338453 1 1 0)
run0042=(0.48847892320123 0 1.0 0.0678338453 1 -1 0)

run_list=(
	run0031
)

#params_file=params_ratio${ratio}.txt
params_file=params_Newtonian_binary.txt

G="1e-12" # 10^{-13}
#G=0

r_inner=2
L=512
N1=64
num="001500"

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

	# text_number=$(printf "%04d" ${run_number})
        subdir=${run}_M${M}_d${d}_mu${mu}_dt_mult${dt_mult}_l${l}_m${m}_Al${Al}_L${L}_N${N1}_complex_rin_${r_inner}_max_level9

	new_dir=Newtonian_${run}_G${G}_max_level9_n${num}_wslope0.25_wradius50
	#new_dir=Newtonian_${run}_G${G}_max_level9_n${num}_noICS
	new_dir_path=${data_directory}/${new_dir}
	#wslope0.25_wradius50
	mkdir -p ${new_dir_path}
	echo "made ${new_dir_path}"

	cd ${work_dir}
	cp pbs_submit_dial3 ${new_dir_path}/pbs_submit
	
	cp ${params_file} ${new_dir_path}/params.txt

	sed -i "s|SUBDIR|${subdir}|" ${new_dir_path}/params.txt
	sed -i "s|GVAL|${G}|" ${new_dir_path}/params.txt
	sed -i "s|OMEGAB|${omega_BH}|" ${new_dir_path}/params.txt
	sed -i "s|OUTNAME|${new_dir}|" ${new_dir_path}/params.txt
	sed -i "s|MUVAL|${mu}|" ${new_dir_path}/params.txt
	sed -i "s|NPLOTFL|${num}|" ${new_dir_path}/params.txt

	cd ${new_dir_path}
	mkdir -p outputs
        cd outputs

	sbatch ../pbs_submit
	#
	cd ${work_dir}
done
