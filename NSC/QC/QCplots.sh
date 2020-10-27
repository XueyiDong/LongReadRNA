#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL


module load R

Rscript QCplots.R