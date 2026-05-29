#!/bin/bash
#SBATCH --job-name=KP10KL0
#SBATCH --partition=astro2_short
#SBATCH --output=KP10KL0.out
#SBATCH --error=KP10KL0.err
#SBATCH --mem=8G

./ANA -analysis "SF" -base_dir "/lustre/astro/semygind/Data/MembraneInclusion/KPKL/KP10KL0/" -init_f 100 -final_f 10000 -final_rep 10