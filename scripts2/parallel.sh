#!/bin/bash
#################### Set the name of the job
#SBATCH --job-name calc_1_1_cauchy
# Launch an array of 100 jobs
## SBATCH --array 1-10
# Specify a time limit
# SBATCH --time 96:00:00
# Redirect stderr and stdout to the same file:
# %A will be replaced by the job ID and %a by the array index
####################
#SBATCH -o calc_1_1_cauchy.out
#SBATCH -e calc_1_1_cauchy.out
# Send email notifications
## SBATCH --mail-type=ALL
# We request an exclusive node for every job in the array
## SBATCH --exclusive
# reserve MB of memory
#SBATCH --mem=60000
# Specify the number of tasks (processes)
#SBATCH --ntasks 1
# Our job is multithreaded, so we ask 4 CPUs per process
#SBATCH --cpus-per-task=31
# So, in total we will have 4x8 CPUs for us
## module load nest/2.10.0
##################################################################
# from here on we can run whatever command we want (e.g., srun python simulate.py $SLURM_ARRAY_JOB_ID)
# we use slurm's environment variables to create unique output files and echo the name of the executing node in that file
##module load anaconda
##source activate r
srun Rscript -e 'source("mccalc_1sym_1snp.R")'
