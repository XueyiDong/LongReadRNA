#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=400G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL


module load R

Rscript prepare_QCplots.R
Rscript QCplots.R
