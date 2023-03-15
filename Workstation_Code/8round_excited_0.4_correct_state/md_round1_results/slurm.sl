#!/bin/bash -l
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=5g
#SBATCH --tmp=5g
#SBATCH --mail-type=NONE
#SBATCH --job-name=ex_MD_1
#SBATCH --mail-user=johan406@umn.edu

module load molpro/2019.2.3
export PATH='/home/goodpast/johan406/iilib/anaconda3/bin':$PATH

stdbuf -oL python < runner.py
