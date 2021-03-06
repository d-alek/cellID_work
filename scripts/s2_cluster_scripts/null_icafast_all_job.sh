#!/bin/bash -l

#########################################
# null ica job: all comps; 10 null mats #
#########################################

#SBATCH --job-name=null_icafast_all_job

#SBATCH --account=punim0613

#SBATCH -p mig

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=1-06:00:00
#SBATCH --mem=100000

#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

#SBATCH --mail-user=adakic@unimelb.edu.au
#SBATCH --mail-type=ALL

# Load required modules and run code
module load R/3.5.2-spartan_gcc-6.2.0

echo "SLURM_JOBID: " $SLURM_JOBID

/usr/bin/time -v Rscript null_icafast_all.r

exit 0

# this resources ok with ~ 6500 x 18500 input data - change accordingly
# .err file will contain max memory required, cpu % usage, run time etc. for re-adjustmets
# .out file will contain n_comp completed for checking progress

# submit job from /punim0613/icafast where the functions are
# data is in /punim0613/icafast/data
# results go to /punim0613/icafast/output