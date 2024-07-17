#!/bin/bash
#SBATCH -a 540,542,544,546,548,550,552,554,556,558,560,562,564,566,568,570,572,574,576,578,580,582,584,586,588,590,592,594,596,598,600
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

julia make_data.jl $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT