#!/bin/bash                
# -- run this code sequentially, srun will give u a node... you could also maybe make this into an sbatch script so u don't have to do it manually but oh well                    

#SBATCH --job-name=build_sys_image
#SBATCH --output=/home/jbenjami/Research_Schneider/jobs/out/sys_img_build_%a_output__%j.out
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --mail-user=jbenjami@caltech.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-node=1


# Get the local path to this script running as script
if [ -n $SLURM_JOB_ID ];  then # `if [ -n $SLURM_JOB_ID ]` checks if $SLURM_JOB_ID is not an empty string
    this_filepath=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')     # check the original location through scontrol and $SLURM_JOB_ID
else
    this_filepath=$(realpath $0)     # otherwise: started with bash. Get the real location.
fi
output=/home/jbenjami/Research_Schneider/jobs/out/sys_img_build_%a_output__%j.out
this_dir=$(dirname ${this_filepath}) # getting location of software_name 
this_file=$(basename $this_dir) # separating the software_name from path

#For running on the node directly
srun --reservation=clima --pty -t 60:00 -N 1 -n 1 --gres=gpu:1 /bin/bash -l
this_filepath=$(pwd)
output=/home/jbenjami/Research_Schneider/jobs/out/sys_img_build.out
this_dir=${this_filepath} # getting location of software_name 
this_file=${this_filepath}/make_image.sh # separating the software_name from path



# echo $this_filepath
# echo $this_dir
# echo $this_file

unset SLURM_STEP_ID
CLIMA=${this_dir}
cd ${CLIMA}
module load julia/1.5.2
module load hdf5/1.10.1 netcdf-c/4.6.1 cuda/10.2 openmpi/4.0.4_cuda-10.2 cmake/3.10.2 # CUDA-aware MPI
julia --project -e 'using Pkg; Pkg.build("MPI")'
echo 'building'
julia --project ${CLIMA}/.dev/systemimage/climate_machine_image.jl ${CLIMA}/ClimateMachine.so
echo 'cleaning up'


mv ${output}  ${this_dir}/sys_img_build_output.out 



