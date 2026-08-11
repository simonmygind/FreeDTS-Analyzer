#!/bin/bash
#SBATCH --job-name=simulation
#SBATCH --partition=astro2_short
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=8G


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#./ANA -analysis SF -base_dir "pathdata" -final_f 10000 -init_f 150 -final_rep 15 -ftmethod DFT -N_k 10
./ANA -analysis SF -base_dir "pathdata" -final_f 4000 -init_f 150 -final_rep 5 7 9 11 13 15 -ftmethod DFT -N_k 10
#./ANA -analysis SF -base_dir "pathdata" -final_f 2000 4000 6000 8000 10000 -init_f 150 -final_rep 15 -ftmethod DFT -N_k 10
echo "Finished simulation
