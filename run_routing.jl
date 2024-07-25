
using JuMP, Gurobi, Plots, Random, CSV, DataFrames, Statistics, Dates, HDF5, LinearAlgebra, FileIO, JLD2, NearestNeighbors

#-----------------------------------------------------------------------------------#

time()

include("scripts/routing/saveandloadassignments.jl")
include("scripts/routing/oldcongestionfunctions.jl")
include("scripts/routing/createroutingnetwork.jl")
include("scripts/routing/solveroutingproblem.jl")

println("Scripts imported")

#const GRB_ENV = Gurobi.Env()

#-----------------------------------------------------------------------------------#

starttime = time()

# Select the instancecd
run_id += 10000 #ifelse(length(ARGS) > 0, parse(Int, ARGS[1]), 1)
#instance_id = 516
#methodname = "dynamicml"
#instancefilename = string(outputfolder, "assignments.jld2") #string("assignments/run", instance_id, "_instance.jld2")
assignmentfilename = string(outputfolder, "assignments.jld2") #string("assignments/run", run_id, "_", methodname, "_assignments.jld2")
outputfilename = string(outputfolder, "routes.csv") #string("outputs/routes_run", run_id, ".csv")

tstep_cong = 15
tstep_r = 30
maxtraveltime = tstep_r * ceil(((warehouse_x_length_meters + warehouse_y_length_meters) / podspeed) / tstep_r)
horizon_r = horizon*2

#-----------------------------------------------------------------------------------#

#Load instance
routenodelookup, routenodes, routenodes_end, numnodes_r, extendednumnodes_r = createroutingnodes()
routearcs, routearclookup, numarcs_r, RA_plus, RA_minus, newarclength, routearcs_path, routearclookup_path, A_space = createroutingarcs()
orders_r, pods_r, podswith_r, workstationassignment, ordersassignedto, orderassignmentbegan, orderopentime, orderclosetime, podsfor, poddelaypenalty, itemsfrom = savetaskassignments(assignmentfilename)
podstartnode_r = routingpodinfo()
podarcset, A_minus_p, A_plus_p, podnodeset = podarcsets_routing(pods_r, newarclength)
congcontribarcs, passingintersections, passingintersectiontimes = getpathcongestioncontributions()
println("Built instance")

#-----------------------------------------------------------------------------------#

#Load assignments
#rders_r, pods_r, podswith_r, workstationassignment, ordersassignedto, orderassignmentbegan = loadtaskassignments(assignmentfilename)
orders_r, pods_r, podswith_r, workstationassignment, ordersassignedto, orderassignmentbegan, orderopentime, orderopentime, podsfor, poddelaypenalty, itemsfrom = savetaskassignments(assignmentfilename)
println("Loaded assignments")

#-----------------------------------------------------------------------------------#

#Solve routing problem
#obj_r, solvetime_r, h_r, y_r, v_r = solveroutingproblem(orders_r, pods_r, podswith_r, workstationassignment, ordersassignedto, orderassignmentbegan)
obj_r, solvetime_r, y_r = solvedelayproblem(orders_r, pods_r, podswith_r, workstationassignment, ordersassignedto, orderassignmentbegan)

#-----------------------------------------------------------------------------------#

#Report out
currpartition = partitioninfo[1]
currsol = partitionsolution[1]
writeroutingmetrics(outputfilename, currsol.itempodpicklist, y_r)

#-----------------------------------------------------------------------------------#

#Visualize
#include("scripts/visualizations/congestionviz.jl")
#congestionviz("nocong.png", 50*warehouse_x_length_meters, 50*(warehouse_y_length_meters), warehouse_x_length_meters, warehouse_y_length_meters, y)

#-----------------------------------------------------------------------------------#

println("Done!")


