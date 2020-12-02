#!/bin/bash                                                                                                                    

#SBATCH --job-name=cfsites_test
#SBATCH --output=/home/jbenjami/Research_Schneider/jobs/%j_%a.out
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:20:00
#SBATCH --mem=0
#SBATCH --mail-user=jbenjami@caltech.edu
#SBATCH --mail-type=ALL

#SBATCH --array=16

set -euo pipefail 

module load julia/1.4.2
module load openmpi/4.0.3_cuda-10.0 cuda/10.0 hdf5/1.10.1 netcdf-c/4.6.1 cmake/3.10.2

export JULIA_MPI_BINARY=system

CLIMA="/home/jbenjami/Research_Schneider/CliMa/ClimateMachine_gcmforcing-cfsite.jl"
julia --project=${CLIMA} -e 'using Pkg; Pkg.instantiate(); Pkg.API.precompile()'

mpirun julia --project=/${CLIMA} ${CLIMA}/experiments/AtmosLES/cfsite_hadgem2-a_07_amip.jl --diagnostics 5smins --output-dir={CLIMA}/GCP_Round2/Site${SLURM_ARRAY_TASK_ID} --group-id=site${SLURM_ARRAY_TASK_ID} --vtk 30smins
