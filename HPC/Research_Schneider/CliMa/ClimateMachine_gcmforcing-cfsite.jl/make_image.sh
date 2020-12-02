#!/bin/bash                
# -- run this code sequentially, srun will give u a node... you could also maybe make this into an sbatch script so u don't have to do it manually but oh well                    


srun --reservation=clima --pty -t 60:00 -N 1 -n 1 --gres=gpu:1 /bin/bash -l
unset SLURM_STEP_ID
CLIMA=/home/jbenjami/Research_Schneider/CliMa/ClimateMachine_gcmforcing-cfsite.jl
cd ${CLIMA}
module load julia/1.5.2
module load hdf5/1.10.1 netcdf-c/4.6.1 cuda/10.0 openmpi/4.0.3_cuda-10.0 cmake/3.10.2 # CUDA-aware MPI
julia --project -e 'using Pkg; Pkg.build("MPI")'
julia --project ${CLIMA}/.dev/systemimage/climate_machine_image.jl ${CLIMA}/ClimateMachine.so