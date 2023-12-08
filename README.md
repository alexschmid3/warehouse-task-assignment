# Amazon Task Assignment - Year 1

This repository implements large-scale neighborhood search (LSNS) for the task
design and scheduling problem with congestion and restricted workload (TDS-CW) in Julia with Gurobi. 

### Basic Use
To solve instances using learn-then-optimize or one of the LSNS benchmarks in the paper, select an or create an instance from `data/indiv_instance_params.csv` and set LSNS parameters in `data/indiv_lsns_parameters.csv`. Set `row_id` to the desired row from the LSNS params file in `run_partitions.jl` and run. 

### Parameter overviews

The size of the warehouse and the number of workstations, pods, items and orders can be set in `data/warehouse_sizes_and_capacities.csv`.  Key control parameters of the LSNS algorithm are as follows:
- `method` - takes values `"LTO"` (learn-then-optimize), `"synergy"` (domain-based heuristic), `"learnbench"` (learning-enhanced benchmark), and `"random"` (randomized)
- `initialization` - takes values `"none"` or `"greedy"`
- `targetnumpods` - target number of pods per LSNS subproblem
- `targetnumorders` - target number of orders per LSNS subproblem
- `targetnumitems` - target number of items per LSNS subproblem
- `subproblemsevaluated`  - number of iterations of LSNS performed
- `subproblemtimelength` - number of time-periods re-optimized per LSNS subproblem
- `timeforreooptimization` - time limit for subproblem re-optimized at each LSNS iteration
- `timeforsubproblemselection`- time limit for subproblem selection at each LSNS iteration
- `tabutype` - time limit for subproblem re-optimized at each LSNS iteration
