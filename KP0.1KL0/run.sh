#!/bin/bash
#SBATCH --job-name=KP0.1KL0
#SBATCH --partition=astro2_short
#SBATCH --output=KP0.1KL0.out
#SBATCH --error=KP0.1KL0.err
#SBATCH --mem=8G

./ANA -analysis "SF" -base_dir "/lustre/astro/semygind/Data/MembraneInclusion/KPKL/KP0.1KL0/" -init_f 100 -final_f 10000 -final_rep 10