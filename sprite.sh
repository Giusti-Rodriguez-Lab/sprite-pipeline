#!/bin/bash
#SBATCH --job-name=sprite
#SBATCH -c 12
#SBATCH -o snakemake.%j.out
#SBATCH -e snakemake.%j.err
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --mem=100G 
#SBATCH --time=100:00:00

module load java
module load snakemake/9.4.0 
module load parabricks/4.5.1 
module load conda/25.1.1

# 2to3 -w misc/sprite-pipeline/python
ml samtools


snakemake -j 1 --cores 12 --use-conda --conda-frontend conda --snakefile Snakefile 
