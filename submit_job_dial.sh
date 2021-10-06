#!/bin/bash
#
# this copy is for the Skylake nodes

work_dir=/home/dc-bamb1/InitialConditionSolver/

output_dir=/scratch/dp092/dc-bamb1/Initial_conditions_solver

L=512
N1=64
box_size=16

data_directory=/p/project/pra116/bamber1/BinaryBHScalarField

# specify the input params for each run I want to submit
# list for each is: mu, delay, dt, G, BH mass ratio

# list for each is: mu, delay, dt, G, BH mass ratio, l, m, Al

run0036=(0.5 0 0.25 0 1 0 0 0)
run0018=(0.5 0 0.25 0.000001 1 0 0 0)
run0023=(0.5 0 0.25 0.0000000001 1 0 0 0) # G = 10^{-10}
run0024=(0.5 0 0.25 0.00000000000000000001 1 0 0 0) # G = 10^{-20}, restart as in run0020
run0026=(0.5 0 0.25 0.000000000000001 1 0 0 0) # G = 10^{-15}
run0027=(0.5 0 0.25 0.00000000000000000001 1 0 0 0) # G = 10^{-20}
run0028=(0.5 0 0.25 0.00000001 1 0 0 0) # G = 10^{-8}
run0029=(1 0 0.0625 0.0000000001 1 0 0 0) # G = 10^{-10}
run0033=(0.5 0 0.25 0.000000001 1 0 0 0) # G = 10^{-9}
run0034=(0.5 0 0.25 0.000000000001 1 0 0 0) # G = 10^{-12}
run0035=(0.5 0 0.25 0.00000000000001 1 0 0 0) # G = 10^{-14}
run0035ICS=(0.5 0 0.25 0.00000000000001 1 0 0 0) # G = 10^{-14}

run_list=(
    run0023
)

for run in "${run_list[@]}"
do
  	cd $work_dir
	# extract parameters
        
        val="$run[0]"; mu="${!val}"
        val="$run[1]"; delay="${!val}"
        val="$run[2]"; dt_mult="${!val}"
        val="$run[3]"; G="${!val}"
        val="$run[4]"; ratio="${!val}"
	 val="$run[5]"; l="${!val}"
        val="$run[6]"; m="${!val}"
        val="$run[7]"; Al="${!val}"

        #if (( $l == 0 ))
        #then
        #   subdir=${run}_mu${mu}_delay${delay}_G${G}_ratio${ratio}
        #else
        #   subdir=${run}_mu${mu}_delay${delay}_G${G}_ratio${ratio}_l${l}_m${m}_Al${Al}
        #fi

	max_level=9
	
	name=${run}_alternative_W_i_decomp_from_Newtonian_file_complex
	#name=${subdir}_initial_conditions
	
        echo ${name} "initial conditions"
        new_dir_path=${output_dir}/${name}
        #
	mkdir -p ${new_dir_path}
        
       	cp pbs_submit_dial ${new_dir_path}/pbs_submit
	params_file=params.txt
        cp ${params_file} ${new_dir_path}/params.txt
        
       	cd ${new_dir_path}
	pwd
        # add the location of the new directory to the params file
        sed -i "s|JOBNAME|${run}IC|" pbs_submit
	sed -i "s|MXLEVEL|${max_level}|" params.txt
	sed -i "s|RUNNAME|${run}_${max_level}_from_Newtonian_file|" params.txt
	sed -i "s|NEWT_FILE|${input_file}|" params.txt
	sed -i "s|GVALUE|${G}|" params.txt
	sed -i "s|MUVALUE|${mu}|" params.txt
	sed -i "s|LSPACE|${L}|" params.txt
	sed -i "s|LSPACE|${L}|" params.txt
        sed -i "s|LSPACE|${L}|" params.txt
        sed -i "s|NBASIC|${N1}|" params.txt
	sed -i "s|NSPACE3|$(($N1/2))|" params.txt
        sed -i "s|BOXSIZE|${box_size}|" params.txt
        sed -i "s|CENTERX|$(($L/2))|" params.txt
        sed -i "s|CENTERY|$(($L/2))|" params.txt
        sed -i "s|CENTERZ|0|" params.txt

        qsub pbs_submit
        #
	cd ${work_dir}
done    


