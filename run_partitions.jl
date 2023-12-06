
#include("run_decomp_dynamicml.jl")

using JuMP, Gurobi, Plots, Random, CSV, DataFrames, Statistics, Dates, HDF5, LinearAlgebra, FileIO, JLD2, NearestNeighbors, SparseArrays

#-----------------------------------------------------------------------------------#

time()
include("scripts/instancegeneration/readdata.jl")
include("scripts/instancegeneration/datageneration.jl")
include("scripts/instancegeneration/initializecongestion.jl")
include("scripts/helper/remove.jl")
include("scripts/helper/sortarcschronologically.jl")
include("scripts/helper/getcurrentobjective.jl")
include("scripts/benchmarks/synergylsns.jl")
include("scripts/benchmarks/findgreedysolution.jl")
include("scripts/decomposition/decomposeproblem.jl")
include("scripts/decomposition/configurepartitioncongestion.jl")
include("scripts/decomposition/writeglobalsolutionoutputs.jl")
include("scripts/lsns/parsetabucodes.jl")
include("scripts/lsns/getmlfeatures.jl")
include("scripts/lsns/enumeratesubproblemwindows.jl")
include("scripts/lsns/createemptysolution.jl")
include("scripts/lsns/updateglobalsolution.jl")
include("scripts/lsns/updatesolvemetrics.jl")
include("scripts/lsns/subproblemselection/selectrandomsubproblem.jl")
include("scripts/lsns/subproblemselection/selectsynergisticsubproblem.jl")
include("scripts/lsns/subproblemselection/selectlearnthenoptimizesubproblem.jl")
include("scripts/lsns/constructsubproblem.jl")
include("scripts/lsns/createemptysolution.jl")
include("scripts/lsns/reoptimizesubproblem.jl")
include("scripts/lsns/updatesolution.jl")
include("scripts/test/checksolution.jl")

println("Scripts imported")

#-----------------------------------------------------------------------------------#

#Debugging mode
visualizationflag = 0
debugprintstatements = 0 			# When a bug is found, set this to 1 for next run to get more detailed error info (warning: it's a lot of print statements, one for each unit test)
debugmode = 0					# 1 --> will perform solution consistency unit tests at each LSNS iteration, use for debugging when changing the code, 0 --> no solution checks, use for computational experiments

#Initialize Gurobi
const GRB_ENV = Gurobi.Env()

# Select the instancecd
row_id = 1 #ifelse(length(ARGS) > 0, parse(Int, ARGS[1]), 1)
instanceparamsfilename = "data/warehouse_sizes_and_capacities.csv"
testingparamsfilename = "data/decomp_instance_parameters.csv"
methodparamsfilename = "data/decomp_lsns_parameters.csv"
instanceparms = CSV.read(instanceparamsfilename, DataFrame)
testingparms = CSV.read(testingparamsfilename, DataFrame)
methodparms = CSV.read(methodparamsfilename, DataFrame)

# Parameter Descriptions:
# ==========================================
# Row 1 = run_id
# Row 2 = test_instance_id
# Row 3 = method
# Row 4 = initialization
# Row 5 = targetnumpods
# Row 6 = targetnumorders
# Row 7 = targetnumitems
# Row 8 = subproblemsevaluated
# Row 9 = subproblemtimeinterval
# Row 10 = timeforreooptimization
# Row 11 = timeforsubproblemselection
# Row 12 = tabutype
# Row 13 = mlmodel

#Get ML training parameters
run_id = methodparms[row_id, 1]
instance_id = methodparms[row_id, 2]
methodname = methodparms[row_id, 3]
solutioninitialization = methodparms[row_id, 4]
targetnumpods = methodparms[row_id, 5]
targetnumorders = methodparms[row_id, 6]
targetnumitems = methodparms[row_id, 7]
numsubproblemsevaluated = 100 #methodparms[row_id, 8]
if methodname == "LTO"
	global numsubproblemsevaluated -= 50
end
subproblemtimelength = methodparms[row_id, 9]
timeforreooptimization = methodparms[row_id, 10]
timeforsubproblemselection = methodparms[row_id, 11]
tabutype = methodparms[row_id, 12]
mlmodelname = methodparms[row_id, 13]
partitionobjective = methodparms[row_id, 14]
minnumpods = targetnumpods - 10
minnumorders = targetnumorders - 10
minnumitems = targetnumitems - 10
numlocalintersections = 6
numhyperlocalintersections = 3 
maxworkstationspersubproblem = 2
windowforcingflag, maxtabu, lastoptpenaltyflag = parsetabucodes(tabutype)

# Parameter Descriptions:
# ==========================================
# Row 1 = instance_id
# Row 2 = warehouse_id (which warehouse layout we're using)
# Row 3 = random seed
# Row 4 = start of time horizon
# Row 5 = end of time horizon

#Get ML training parameters from file
#instance_id = testingparms[run_id, 1]
warehouse_id = testingparms[instance_id, 2]
random_seed = testingparms[instance_id, 3]
horizonstart = testingparms[instance_id, 4]
horizonend = testingparms[instance_id, 5]
horizon = horizonend - horizonstart
stationsperpartition = testingparms[instance_id, 6]

# Parameter Descriptions:
# ==========================================
# Row 1 = warehouse_id
# Row 2 = warehouse_x_length_meters
# Row 3 = warehouse_y_length_meters
# Row 4 = workstations
# Row 5 = unique_items
# Row 6 = items_per_pod
# Row 7 = orders_per_hour
# Row 8 = prob_one_item_orders
# Row 9 = geo_dist_param_orders
# Row 10 = time_step
# Row 11 = congestion_time_step
# Row 12 = pod_process_time (# of seconds it takes to process a new pod at the workstation)
# Row 13 = item_process_time (# of seconds it takes to pick an item at the workstation))
# Row 14 = workstation_order_capacity
# Row 15 = intersection_capacity

#Get warehouse parameters from file
warehouse_x_length_meters = instanceparms[warehouse_id, 2]
warehouse_y_length_meters = instanceparms[warehouse_id, 3]
num_workstations = instanceparms[warehouse_id, 4]
num_unique_items = instanceparms[warehouse_id, 5]
num_items_per_pod = instanceparms[warehouse_id, 6]
num_orders_per_hour = instanceparms[warehouse_id, 7]
prob_one_item_orders = instanceparms[warehouse_id, 8]
geo_dist_param_orders = instanceparms[warehouse_id, 9]
tstep = instanceparms[warehouse_id, 10]
congestiontstep = instanceparms[warehouse_id, 11]
podprocesstime = instanceparms[warehouse_id, 12]
itemprocesstime = instanceparms[warehouse_id, 13]
workstationordercapacity = instanceparms[warehouse_id, 14]
intersectioncapacity = instanceparms[warehouse_id, 15]
capacitybuffer = 10
podspeed = 1 #Current assumption: 1 meter per second, ~2.2 miles per hour
generation_warmstart_flag = 0

println("Parameters read")

#Files
mlmodelfilename = string("models/", mlmodelname, ".jld2")
outputfolder = string("outputs/run", run_id,"_", today())
globalsolutionfilename = string(outputfolder, "/output.csv")
if !(isdir(outputfolder))
	mkdir(outputfolder)
end
visualizationfolder = string(outputfolder, "/viz")
if !(isdir(visualizationfolder))
	mkdir(visualizationfolder)
end

#Initialize timer
time()

#-----------------------------------------------------------------------------------#

#Create the instance
println("----------------------Create instance---------------------")

#Call data generation function
instancedetails = Data(warehouse_x_length_meters, warehouse_y_length_meters, num_workstations, num_unique_items, num_items_per_pod, num_orders_per_hour, prob_one_item_orders, geo_dist_param_orders, random_seed)
vertstreets, horstreets, workloc, podstart = instancedetails[1][1], instancedetails[1][2], instancedetails[1][3], instancedetails[1][4]
inventoryinfo = instancedetails[2]
orderheader, orderdetail = instancedetails[3][1], instancedetails[3][2]

println("Instance randomized")

#Set up problem and relevant parameters for input to optimization model
storagelocs, workstations, loccoords, loclookup, numstoragelocs, numworkstations = getlocations(podstart, workloc)
intersections, queueintersections, intcoords, intlookup, intcoords_nn = getintersections(horstreets, vertstreets, numstoragelocs, numworkstations)
maploctointersection = maplocationaccesspoints(storagelocs, workstations, intcoords_nn, loccoords, intlookup)
numnodes, nodes, nodelookup, N_end, times = createtimespacenetwork(storagelocs, workstations, [])
extendednumnodes, extendednodes, nodelookup, N_ext_before, N_ext_after, extendedtimes, dummystarttime, dummyendtime = extendtimespacenetwork(nodelookup)
numpods, allpods, podstartnode, podstorageloc = getpods(podstart)
podswith, allitems, inventory, podstartinventory = getinventory(inventoryinfo)
orders, itemson, deadline, ordersdue = getorders(orderheader, orderdetail, allitems, horizonstart, horizonend)
podsstoredat = getpodsstoredat()
C = getobjective()
warehousedistance = calcwarehousedistances()

println("Instance translated")

prearcs, arclength = getphysicalarcs(podspeed, loccoords)
numarcs, arcs, arclookup, A_plus, A_minus, A_space, A_queues = createtimespacearcs(prearcs, numnodes) #Slow
extendednumarcs, extendedarcs, arclookup, A_plus, A_minus, A_space = extendtimespacearcs(arclookup, A_space, A_plus, A_minus)

#Changes orders with multiple of the same item - temp fix, though shouldn't impact results
itemson, items, orderswith = orderitemtempfix(itemson)

#Needed for Riley's outputs: Finish problem set up once the order items have been updated
pods, podarcset, A_minus_p, A_plus_p, podnodeset, arcpodset = podarcsets(items, arclength) #Slower

println("Instance prepared for optimization")

#-----------------------------------------------------------------------------------#

if visualizationflag == 1
	include("scripts/visualizations/warehouseviz.jl")
	include("scripts/visualizations/workstationviz.jl")
end

#Helper function
function arcDesc(a)
	println(nodelookup[arclookup[a][1]], " ==> ", nodelookup[arclookup[a][2]])
end

#-----------------------------------------------------------------------------------#

#Congestion initialization
currcong, maps, congestionsignature, intersectionmaxpods, intersectiontimemaxpods = initializecongestion()

#-----------------------------------------------------------------------------------#

#Read the desired ML model and get problem features
beta, features, featuresfor, featurenums, featureinfo = getmlfeatures(mlmodelfilename)

#-----------------------------------------------------------------------------------#

println("--------------------Decompose instance-------------------")

#Warehouse visualization and statistics
if visualizationflag == 1
	warehouseviz(string(outputfolder, "/warehouselayout.png"), 4000)
	orderlevels = [length([m for m in orders if length(itemson[m]) == l]) for l in 1:20]
	plot1 = Plots.bar(1:20, orderlevels)
	savefig(plot1,string(outputfolder, "/ordersize.png"))

	invlevels = [length([i for i in items if length(podswith[i]) == l]) for l in 1:25]
	plot2 = Plots.bar(1:25, invlevels)
	savefig(plot2,string(outputfolder, "/poditemdistribution.png"))
end

#Decompose problem into partitions
numpartitions, partitions, partitioninfo, globalpartitionid = decomposeproblem(stationsperpartition, partitionobjective, beta, features, featureinfo, featurenums)
globalpartition = partitioninfo[globalpartitionid]

#Enumerate windows for the global partition
windows_global, windowsduring_global, windowidlookup_global, windowscontaining_global = enumeratesubproblemwindows(globalpartition, 1, subproblemtimelength)
windowsynergy_global = preprocesswindowsynergies(windows_global, windowidlookup_global, featureinfo, features, beta, featurenums)

#-----------------------------------------------------------------------------------#

println("----------------------Solve instance----------------------")

#Initialize global solution
globalsolution = createemptysolution(orders, pods, podswith, workstations)
partitionsolution = Dict()
for s in 1:numpartitions
	currpartition = partitioninfo[s]
	partitionsolution[s] = createemptysolution(currpartition.orders, currpartition.pods, currpartition.podswith, currpartition.workstations)
end
solvemetrics = (solve_time=zeros(numpartitions+1), solvetime_init=zeros(numpartitions+1), solvetime_spsel=zeros(numpartitions+1), solvetime_sp=zeros(numpartitions+1), lsnsiterations=zeros(numpartitions+1))
writeglobalsolutionoutputs_init(globalsolutionfilename)
counter = 1

#Solve each partition
for s in 1:numpartitions

	println("===== PARTITION $s =====")

	#Get partition info
	currpartition = partitioninfo[s]

	#Find an initial solution
	currsol = createemptysolution(currpartition.orders, currpartition.pods, currpartition.podswith, currpartition.workstations)
	initstarttime = time()
	if solutioninitialization == "greedy"
		currsol = findgreedysolution(currpartition, currsol, currpartition.orders, currpartition.pods, currpartition.podswith, currpartition.workstations)
	elseif solutioninitialization == "biggreedy"
		currsol = findgreedysolution(currpartition, currsol, [m for m in currpartition.orders if length(itemson[m]) >= 4], currpartition.pods, currpartition.podswith, currpartition.workstations)
	end
	initializationtime = time() - initstarttime

	#Report on initial solution
	getcurrentobjective(currpartition, currsol)
	if visualizationflag == 1
		workstationviz(string(visualizationfolder,"/station_partition", s,"_initial.png"), currpartition, currsol)
	end
	writeglobalsolutionoutputs_iter("init", initializationtime, initializationtime, 0, 0, globalsolutionfilename, s, currpartition, currsol)
	solvemetrics.solvetime_init[s] += initializationtime

	#Find subproblem windows
	windows, windowsduring, windowidlookup, windowscontaining = enumeratesubproblemwindows(currpartition, maxworkstationspersubproblem, subproblemtimelength)
	windowsynergy = preprocesswindowsynergies(windows, windowidlookup, featureinfo, features, beta, featurenums)

	#Initialize LSNS sets
	inittime = time()
	tabulist, lastoptimizeddifference = [], zeros(length(windows))

	#Improve solution iteratively via LSNS
	for sp_iter in 1:numsubproblemsevaluated

		println("----- ITER $sp_iter -----")
		iterationstarttime = time()

		#Select subproblem elements
		spselectionstarttime = time()
		if methodname == "LTO"
			sp_winid, sp_orders, sp_window, sp_pods, sp_itemson, sp_items, tabulist, predicted_obj = selectlearnthenoptimizesubproblem(currpartition, currsol, windows, windowscontaining, windowidlookup, windowsduring, windowsynergy, targetnumorders, targetnumpods, tabulist, lastoptimizeddifference, sp_iter)
		elseif methodname == "random"
			sp_winid, sp_orders, sp_window, sp_pods, sp_itemson, sp_items = selectrandomsubproblem(currpartition, windows, windowidlookup, currsol, targetnumorders, targetnumpods)
			predicted_obj = 0
		elseif methodname == "synergy"
			sp_winid, sp_orders, sp_window, sp_pods, sp_itemson, sp_items, tabulist = selectsynergisticsubproblem(currpartition, windows, currsol, targetnumorders, targetnumpods, rand(1:maxworkstationspersubproblem), tabulist)
			predicted_obj = 0
		end
		spselectiontime = time() - spselectionstarttime
		println("Selection time = ", spselectiontime, " seconds")

		#Build subproblem structure
		spconstructstarttime = time()
		sp = constructsubproblem(currpartition, sp_orders, sp_window, sp_pods, sp_itemson, sp_items, currsol)
		println("Construction time = ", time() - spconstructstarttime, " seconds")

		#Re-optimize subproblem
		sp_obj, sp_solvetime, h_sp, y_sp, z_sp, f_sp, g_sp, v_sp, feasibleflag_sp = reoptimizesubproblem(sp, currsol, currpartition)
		println("Re-opt time = ", sp_solvetime, " seconds")

		#Update solution
		spupdatestarttime = time()
		currsol = updatesolution(sp, currsol, currpartition, sp_obj, h_sp, y_sp, v_sp)
		lastoptimizeddifference = updatelastoptimizeddifference(lastoptimizeddifference, tabulist, sp_winid, windows, predicted_obj, sp_obj)
		iterationtime = time() - iterationstarttime
		println("Update time = ", time() - spupdatestarttime, " seconds")
		updatesolvemetrics(s, iterationtime, sp_solvetime, spselectiontime)
		println("Total time = ", iterationtime, " seconds")

		#Write solution metrics
		writeglobalsolutionoutputs_iter(sp_iter, 0, iterationtime, spselectiontime, sp_solvetime, globalsolutionfilename, s, currpartition, currsol)

		#Visualize new solution
		if visualizationflag == 1
			workstationviz(string(visualizationfolder, "/station_partition", s,"_iter", sp_iter,".png"), currpartition, currsol)
		end

		#For debugging: check consistency of current solution
		if debugmode == 1
			checksolution(currpartition, currsol, debugprintstatements)
		end

		println("Partition throughput = ", sum(length(currsol.itempodpicklist[w,t]) for w in currpartition.workstations, t in times))

		global counter += 1

	end

	#Update global solution
	updateglobalsolution(currpartition, currsol)
	partitionsolution[s] = currsol

end

#-----------------------------------------------------------------------------------#

writeglobalsolutionoutputs(globalsolutionfilename, solvemetrics)
