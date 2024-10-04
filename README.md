# Amazon Year 1 

This repository implements large-scale neighborhood search (LSNS) for the task
design and scheduling problem with congestion and restricted workload (TDS-CW) in Julia with Gurobi. The corresponding paper is ["Robotic warehousing operations: a learn-then-optimize approach to large-scale neighborhood search"](https://arxiv.org/abs/2408.16890).

The paper proposes a novel learn-then-optimize approach to guide subproblem selection in LSNS. 

### Basic Use
Run `run_partitions.jl` to solve an individual instance of the TDS-CW with LSNS (either via learn-then-optimize or another benchmark).  Set `row_id` to the desired row from the LSNS parameters file (`data/test_run_parameters.csv`) in `run_partitions.jl` and run. 

### Data and parameters

All data for the project is synthetic. The size of the warehouse and the number of workstations, pods, items and orders are set in `data/warehouse_sizes_and_capacities.csv`. Parameters of individual instances of the TDS-CW are stored in `data/test_instance_params.csv`. Finally, algorithmic controls of the LSNS algorithm are stored in `data/test_run_parameters.csv`.

 Key control parameters of the LSNS algorithm are as follows:
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
