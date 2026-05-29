#!/bin/bash
#SBATCH --job-name=KP10KL5
#SBATCH --partition=astro2_short
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --mem=8G

./ANA -analysis "SF" -base_dir "/lustre/astro/semygind/Data/MembraneInclusion/KPKL/KP10KL5/" -final_f 10000 -final_rep 10 -init_f 100