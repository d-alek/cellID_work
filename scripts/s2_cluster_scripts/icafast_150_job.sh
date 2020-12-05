#!/bin/bash -l

#########################################
# icafast loop-job: 150 comps; 100 reps #
#########################################

#SBATCH --job-name=icafast_150_job

#SBATCH --account=punim0613

#SBATCH -p snowy

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-06:00:00
#SBATCH --mem=70000

#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

#SBATCH --mail-user=adakic@unimelb.edu.au
#SBATCH --mail-type=ALL

# Load required modules and run code
module load R/3.5.2-spartan_gcc-6.2.0

echo "SLURM_JOBID: " $SLURM_JOBID

/usr/bin/time -v Rscript icafast_150.r 

exit 0

# this resources ok with ~ 6500 x 18500 input data matrix - change accordingly
# .err file will contain max memory required, cpu % usage, run time etc. for re-adjustmets
# .out file will contain n_reps completed for checking progress

# submit job from /punim0613/icafast where the scripts are
# data is in /punim0613/icafast/data
# results go to /punim0613/icafast/output