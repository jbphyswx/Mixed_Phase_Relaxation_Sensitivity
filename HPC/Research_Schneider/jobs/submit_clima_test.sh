#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=00:29:59   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=4G   # memory per CPU core
#SBATCH -J "jb_clima_test"   # job name
#SBATCH --mail-user=jbenjami@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --qos=debug
#SBATCH --reservation=clima
#SBATCH --output=/home/jbenjami/Research_Schneider/jobs/test_rb.out

#set -euo pipefail
#set -x

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load julia/1.4.2
module load cuda/10.0 openmpi/4.0.3_cuda-10.0 hdf5/1.10.1 netcdf-c/4.6.1 

export JULIA_DEPOT_PATH=/home/jbenjami/.julia
export JULIA_MPI_BINARY=system 
export OPENBLAS_NUM_THREADS=1
export PATH="/usr/sbin:$PATH"

julia --project=/home/jbenjami/Research_Schneider/ClimateMachine.jl /home/jbenjami/Research_Schneider/ClimateMachine.jl/tutorials/Atmos/dry_rayleigh_benard.jl
#julia --project /home/jbenjami/Research_Schneider/ClimateMachine.jl/tutorials/Atmos/dry_rayleigh_benard.jl
