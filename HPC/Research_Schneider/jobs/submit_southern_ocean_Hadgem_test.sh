#!/bin/bash                                                                                                                    

#SBATCH --job-name=cfsites_test
#SBATCH --output=/home/jbenjami/Research_Schneider/jobs/out/site_%a_output__%j.out
#SBATCH --nodes=1
#SBATCH --time=05:30:00
#SBATCH --mail-user=jbenjami@caltech.edu
#SBATCH --mail-type=ALL
#SBATCH --array=82
#SBATCH --ntasks=1
#SBATCH --gres=gpu:2
#SBATCH --ntasks-per-node=1

# can also do stuff like #SBATCH --array=17-23
dt=$(date +"%Y_%m_%d__%H%M_%S")
data_path="/home/jbenjami/Research_Schneider/Data/cfsites/CMIP5/CFMIP2/southern_ocean_forcing/"
model="HadGEM2-A"
exper="amip"
rip="r1i1p1"
# years="all"

# # pretty good if you intend to addd some surface heating and midlevel cooling, predates big convective event by a few days
# years=[2006]
# months=[2]
# days=[12]
# hours=[19]
# minutes=[00]

# is in the midle of convection, has no shot at instasbility...
# years=[2006]
# months=[2]
# days=[14]
# hours=[9]
# minutes=[30]

# # very unstable at bottom, kinda unstable up to around 6000m or so
# years=[2006]
# months=[2]
# days=[05]
# hours=[05]
# minutes=[30]

# # Unstable but only up to maybe 6km or so but violently so at surface....
years=[2006]
months=[2]
days=[24]
hours=[2]
minutes=[00]

# # midlevel warm spot at around 5 K/km to overcome
# years=[2006]
# months=[2]
# days=[29]
# hours=[09]
# minutes=[30]

delta_h=1000
delta_v=75
xmax=20000
ymax=20000
zmax=15000 # something like 15000 will be the goal
tmax=$((3600*12))
# tmax=$((3600*12))
# tmax=$((6))
tau_sub_dep_scale=5000
tau_cond_evap_scale=5
solver_type="imex_solver"
# solver_type="LSRK144NiegemannDiehlBusch"
# solver_type="imex_solver"
# timestep=.05
moisture_model="nonequilibrium"
# moisture_model="equilibrium"

# precipitation_model="rain"
precipitation_model="remove_precipitation"
# precipitation_model="no_precipitation"

set -euo pipefail # kill the job if anything fails
# set -x # echo script line by line as it runs it

# # Select version of CliMa you want
# CLIMA=/home/jbenjami/Research_Schneider/CliMa/ClimateMachine_relax.jl
CLIMA=/home/jbenjami/Research_Schneider/CliMa/_ClimateMachine_gcmforcing-cfsite.jl
# CLIMA=/home/jbenjami/Research_Schneider/CliMa/ClimateMachine_gcmforcing-cfsite.jl
# CLIMA=/home/jbenjami/Research_Schneider/CliMa/ClimateMachine.jl

in_script=${CLIMA}/experiments/AtmosLES/cfsite_hadgem2-a_07_amip.jl
outdir=/home/jbenjami/Research_Schneider/Data/cfsites/CMIP5/CFMIP2/output/so_hadgem_test/site_${SLURM_ARRAY_TASK_ID}_${dt}/

# SEE https://github.com/CliMA/ClimateMachine.jl/wiki/Caltech-Central-Cluster#sample-sbatch-scripts

# # Run normally w/ GPU
module purge
# module load julia/1.4.2
# module load openmpi/4.0.3_cuda-10.0 cuda/10.0 hdf5/1.10.1 netcdf-c/4.6.1 cmake/3.10.2
module load julia/1.5.2
module load hdf5/1.10.1 netcdf-c/4.6.1 cuda/10.2 openmpi/4.0.4_cuda-10.2 cmake/3.10.2 # CUDA-aware MPI
# I think these 3 gotta roll as a package
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
export JULIA_DEPOT_PATH=${CLIMA}/.julia_depot # not sure why i commented this out, seems to lead to mpi binary has changed, please run pkg.build mpi...
export JULIA_MPI_BINARY=system # does this cause a problem between nodes?
export JULIA_CUDA_USE_BINARYBUILDER=false


# # # Run on CPU (also set #gpu to 0)
# module purge
# module load julia/1.5.2
# module load hdf5/1.10.1 netcdf-c/4.6.1 openmpi/4.0.1
# export JULIA_MPI_BINARY=system
# export JULIA_CUDA_USE_BINARYBUILDER=false
# export JULIA_DEPOT_PATH=$CLIMA/.julia_depot


# see https://stackoverflow.com/questions/56962129/how-to-get-original-location-of-script-used-for-slurm-job
# check if script is started via SLURM or bash, if with SLURM: there variable '$SLURM_JOB_ID' will exist
if [ -n $SLURM_JOB_ID ];  then # `if [ -n $SLURM_JOB_ID ]` checks if $SLURM_JOB_ID is not an empty string
    this_filepath=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')     # check the original location through scontrol and $SLURM_JOB_ID
else
    this_filepath=$(realpath $0)     # otherwise: started with bash. Get the real location.
fi

# thisdir=$(dirname $(dirname ${this_filepath})) # getting location of software_name , call twice because 
# thisfile=$(basename $thisdir) # separating the software_name from path

# thisdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )" << Doesn't work in slurm, "unbound variable bash_source" error
# thisfile=$(basename ${BASH_SOURCE[0]})
# this_filepath=${thisdir}/${thisfile}
mkdir ${outdir}
cp ${this_filepath} ${outdir}
cp ${in_script}     ${outdir}

# # Akshay way
# julia --project=${CLIMA} -e 'using Pkg; Pkg.instantiate(); Pkg.API.precompile()'

# see https://github.com/CliMA/ClimateMachine.jl/wiki/Caltech-Central-Cluster#sample-sbatch-scripts
julia --project=${CLIMA} -e 'using Pkg; Pkg.instantiate(); Pkg.build()' 
julia --project=${CLIMA} -e 'using Pkg; Pkg.precompile()'

# # can also try to run off a system image a la https://clima.github.io/ClimateMachine.jl/latest/DevDocs/SystemImage/, not sure is faster, also have to be careful wit imports u choose
# julia -J${CLIMA}/ClimateMachine.so  --project # honestly seems no faster


 # think it's supposed to be a dashes here, gets resolved to _ in julia bc we told it to lmao
 # comment things out inside backticks, as per https://stackoverflow.com/questions/9522631/how-to-put-a-line-comment-for-a-multi-line-command or maybe switch to array syntax (same link) later...
args=(--project=/${CLIMA} 
      ${in_script} 
      --diagnostics 5smins 
      --output-dir=${outdir} 
      --data-path=${data_path} 
      --model=${model} 
      --exper=${exper} 
      --rip=${rip} 
      --years=${years} 
      --months=${months} 
      --days=${days} 
      --hours=${hours} 
      --minutes=${minutes} 
      --sites=${SLURM_ARRAY_TASK_ID} 
      --delta-h=${delta_h} 
      --delta-v=${delta_v} 
      --xmax=${xmax} 
      --ymax=${ymax} 
      --zmax=${zmax} 
      --tmax=${tmax} 
      # --timestep=${timestep} 
      --tau-cond-evap-scale=${tau_cond_evap_scale} 
      --tau-sub-dep-scale=${tau_sub_dep_scale} 
      --moisture-model=${moisture_model} 
      --precipitation-model=${precipitation_model}
      --solver-type=${solver_type} 
      --vtk 30smins
      )
# mpirun julia "${args[@]}" 
mpiexec julia "${args[@]}" # allegedly more standardized

# can also try to run off a system image a la https://clima.github.io/ClimateMachine.jl/latest/DevDocs/SystemImage/, not sure is faster, also have to be careful wit imports u choose
# mpiexec julia -J${CLIMA}/ClimateMachine.so "${args[@]}" # --project


# mpirun julia --project=/${CLIMA} \
#   ${in_script} \
#   --diagnostics 5smins \
#   --output-dir=${outdir} \
#   --data-path=${data_path} \
#   --model=${model} \
#   --exper=${exper} \
#   --rip=${rip} \
#   --years=${years} \
#   --months=${months} \
#   --days=${days} \
#   --hours=${hours} \
#   --minutes=${minutes} \
#   --sites=${SLURM_ARRAY_TASK_ID} \
#   --delta-h=${delta_h} \
#   --delta-v=${delta_v} \
#   --xmax=${xmax} \
#   --ymax=${ymax} \
#   --zmax=${zmax} \
#   --tmax=${tmax} \
#   `# --timestep=${timestep}` \ 
#   --tau-cond-evap-scale=${tau_cond_evap_scale} \
#   --tau-sub-dep-scale=${tau_sub_dep_scale} \
#   --moisture-model=${moisture_model} \
#   --solver-type=${solver_type} \
#   --vtk 30smins

# TO DO: add in a try except to still do some cleanup and not leave dangling folders if the run fails

echo '... cleaning up ...'
# # see https://stackoverflow.com/questions/56962129/how-to-get-original-location-of-script-used-for-slurm-job
# # check if script is started via SLURM or bash, if with SLURM: there variable '$SLURM_JOB_ID' will exist
# if [ -n $SLURM_JOB_ID ];  then # `if [ -n $SLURM_JOB_ID ]` checks if $SLURM_JOB_ID is not an empty string
#     this_filepath=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')     # check the original location through scontrol and $SLURM_JOB_ID
# else
#     this_filepath=$(realpath $0)     # otherwise: started with bash. Get the real location.
# fi

# thisdir=$(dirname $(dirname ${this_filepath})) # getting location of software_name 
# thisfile=$(basename $thisdir) # separating the software_name from path

# # thisdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )" << Doesn't work in slurm, "unbound variable bash_source" error
# # thisfile=$(basename ${BASH_SOURCE[0]})
# # this_filepath=${thisdir}/${thisfile}

# cp ${this_filepath} ${outdir}
# cp ${in_script}     ${outdir}
mv /home/jbenjami/Research_Schneider/jobs/out/site_${SLURM_ARRAY_TASK_ID}_output__${SLURM_JOB_ID}.out  ${outdir}/site_${SLURM_ARRAY_TASK_ID}_${dt}_output__${SLURM_JOB_ID}.out 
echo '... Done!'