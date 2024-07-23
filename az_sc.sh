#!/bin/bash
#SBATCH -a 402,404,406,408,410,412,414,416,418,420
#SBATCH --partition=xeon-p8
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=0-06:00
#SBATCH -o /home/gridsan/aschmid/warehouse-task-assignment/outerr/mainrun_%a.out
#SBATCH -e /home/gridsan/aschmid/warehouse-task-assignment/outerr/mainrun_%a.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=aschmid@mit.edu

module load julia/1.9.2
module load gurobi/gurobi-1102

julia run_partitions.jl $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT