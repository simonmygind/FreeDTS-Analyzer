#!/bin/bash
#SBATCH --job-name=KP0KL0
#SBATCH --partition=astro2_short
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --mem=8G

./ANA -analysis "SF" -base_dir "/lustre/astro/semygind/Data/MembraneInclusion/KPKL/KP0KL0/" -final_f 10000 -final_rep 10 -init_f 100