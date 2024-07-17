#!/bin/bash
#SBATCH -a 480,482,484,486,488,490,492,494,496,498,500,502,504,506,508,510,512,514,516,518,520,522,524,526,528,530,532,534,536,538
#SBATCH --partition=xeon-p8
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=0-03:00
#SBATCH -o /home/gridsan/aschmid/warehouse-task-assignment/outerr/maintrain_%a.out
#SBATCH -e /home/gridsan/aschmid/warehouse-task-assignment/outerr/maintrain_%a.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=aschmid@mit.edu

module load julia/1.9.2
module load gurobi/gurobi-1102

julia run_makedata.jl $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT