#!/bin/bash
# Runs 04b_WGCNA_signed_network.R on the Rutgers Amarel cluster.
# Submit from the repository root:  sbatch 04_WGCNA/04b_run_WGCNA_amarel.sh
#SBATCH --job-name=wgcna_trans25
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=wgcna_%j.out
#SBATCH --error=wgcna_%j.err

module load apptainer/1.2.5

# R 4.5.1 container provided by Amarel
SIF=/projects/community/singularity.images/rstudio/R4.5.1_Rstudio_updated.sif

# R packages (WGCNA, DESeq2, dplyr) were installed to this user library
export R_LIBS_USER=/scratch/tab307/Rlibs

apptainer exec --bind /scratch,/home "$SIF" Rscript 04_WGCNA/04b_WGCNA_signed_network.R
