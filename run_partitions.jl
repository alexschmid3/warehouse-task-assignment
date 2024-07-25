
using JuMP, Gurobi, Plots, Random, CSV, DataFrames, Statistics, Dates, HDF5, LinearAlgebra, FileIO, JLD2, NearestNeighbors, SparseArrays

#-----------------------------------------------------------------------------------#

include("scripts/instancegeneration/readdata.jl")
include("scripts/instancegeneration/datageneration.jl")
include("scripts/instancegeneration/initializecongestion.jl")
include("scripts/helper/remove.jl")
include("scripts/helper/sortarcschronologically.jl")
include("scripts/helper/getcurrentobjective.jl")
include("scripts/helper/parsemethodname.jl")
include("scripts/greedy/findgreedysolution.jl")
include("scripts/decomposition/decomposeproblem.jl")
include("scripts/decomposition/configurepartitioncongestion.jl")
include("scripts/decomposition/writeglobalsolutionoutputs.jl")
include("scripts/decomposition/updateglobalsolution.jl")
include("scripts/lsns/parsetabucodes.jl")
include("scripts/lsns/enumeratesubproblemwindows.jl")
include("scripts/lsns/createemptysolution.jl")
include("scripts/lsns/updatesolvemetrics.jl")
include("scripts/lsns/constructsubproblem.jl")
include("scripts/lsns/reoptimizesubproblem.jl")
include("scripts/lsns/updatesolution.jl")
include("scripts/lsns/modelfeatures/getlearningbenchmarkfeatures.jl")
include("scripts/lsns/modelfeatures/preprocesslearningbenchmarkfeatures.jl")
include("scripts/lsns/modelfeatures/getmlfeatures.jl")
include("scripts/lsns/modelfeatures/updatecongestionfeatures.jl")
include("scripts/lsns/subproblemselection/selectrandomsubproblem.jl")
include("scripts/lsns/subproblemselection/selectsynergisticsubproblem.jl")
include("scripts/lsns/subproblemselection/selectlearnthenoptimizesubproblem.jl")
include("scripts/lsns/subproblemselection/selectlearningbenchmarksubproblem.jl")
include("scripts/test/checksolution.jl")
include("scripts/figures/histograms.jl")

#include("scripts/visualizations/warehouseviz.jl")

println("Scripts imported")

#-----------------------------------------------------------------------------------#

#Toggle switches
visualizationflag = 0		       # 1 --> Produce visualizations of warehouse and workstation solutions
debugmode = 0					   # 1 --> will perform solution consistency unit tests at each LSNS iteration, use for debugging when changing the code; 0 --> no solution checks, use for computational experiments
debugprintstatements = 0 		   # When a bug is found, set this to 1 for next run to get more detailed error info (warning: it's a lot of print statements, one for each unit test)
subproblemstatsreporting_flag = 0
ordergraphreporting_flag = 0

#Initialize Gurobi
const GRB_ENV = Gurobi.Env()

# Select the run files
row_id = ifelse(length(ARGS) > 0, parse(Int, ARGS[1]), 1) # (for cluster submissions)
warehouseparamsfilename = "data/warehouse_sizes_and_capacities.csv"
instanceparamsfilename = "data/decomp_instance_parameters.csv"
methodparamsfilename = "data/decomp_run_parameters.csv" #extensions/orderslots/
projectfolder = "outputs/partitionruns/" #"outputs/routing/"
warehouseparms = CSV.read(warehouseparamsfilename, DataFrame)
instanceparms = CSV.read(instanceparamsfilename, DataFrame)
methodparms = CSV.read(methodparamsfilename, DataFrame)

nocongestion_flag = 0 #methodparms[row_id, 18]

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
# Row 9 = subproblemtimelength
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
numsubproblemsevaluated = methodparms[row_id, 8]
subproblemtimelength = methodparms[row_id, 9]
timeforreooptimization = methodparms[row_id, 10]
timeforsubproblemselection = methodparms[row_id, 11]
tabutype = methodparms[row_id, 12]
mlmodelname = methodparms[row_id, 13]
partitionobjective = methodparms[row_id, 14]
stationtostation_flag = methodparms[row_id, 15]
targetnumworkstations = methodparms[row_id, 16]
anystoragelocation_flag = methodparms[row_id, 17]
minnumpods = targetnumpods - 10
minnumorders = targetnumorders - 10
minnumitems = targetnumitems - 10
numlocalintersections = 6
numhyperlocalintersections = 3 
maxworkstationspersubproblem = 2

windowforcingflag, maxtabu, lastoptpenaltyflag = parsetabucodes(tabutype)
shortmethodname, subproblembudget = parsemethodname(methodname)

# Parameter Descriptions:
# ==========================================
# Row 1 = instance_id
# Row 2 = warehouse_id (which warehouse layout we're using)
# Row 3 = random seed
# Row 4 = start of time horizon
# Row 5 = end of time horizon

#Get ML training parameters from file
warehouse_id = instanceparms[instance_id, 2]
random_seed = instanceparms[instance_id, 3]
horizonstart = instanceparms[instance_id, 4]
horizonend = instanceparms[instance_id, 5]
horizon = horizonend - horizonstart
stationsperpartition = instanceparms[instance_id, 6]

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
warehouse_x_length_meters = warehouseparms[warehouse_id, 2]
warehouse_y_length_meters = warehouseparms[warehouse_id, 3]
num_workstations = warehouseparms[warehouse_id, 4]
num_unique_items = warehouseparms[warehouse_id, 5]
num_items_per_pod = warehouseparms[warehouse_id, 6]
num_orders_per_hour = warehouseparms[warehouse_id, 7]
prob_one_item_orders = warehouseparms[warehouse_id, 8]
geo_dist_param_orders = warehouseparms[warehouse_id, 9]
tstep = warehouseparms[warehouse_id, 10]
congestiontstep = warehouseparms[warehouse_id, 11]
podprocesstime = warehouseparms[warehouse_id, 12]
itemprocesstime = warehouseparms[warehouse_id, 13]
workstationordercapacity = warehouseparms[warehouse_id, 14]
intersectioncapacity = warehouseparms[warehouse_id, 15]
capacitybuffer = 10
podspeed = 1 #Current assumption: 1 meter per second, ~2.2 miles per hour
generation_warmstart_flag = 0

println("Parameters read")

#Files
mlmodelfilename = string("models/", mlmodelname, ".jld2")
outputfolder = string(projectfolder,"run", run_id,"_", today())
if !(isdir(projectfolder))
	mkdir(projectfolder)
end
globalsolutionfilename = string(outputfolder, "/output.csv")
if !(isdir(outputfolder))
	mkdir(outputfolder)
end
visualizationfolder = string(outputfolder, "/viz")
if !(isdir(visualizationfolder))
	mkdir(visualizationfolder)
end
subproblemstatsreportingfilename = string(outputfolder, "/subproblems.csv")
ordergraphreportingfilename = string(outputfolder, "/ordergraph.csv")

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
if shortmethodname == "learnbench"
	beta = load(mlmodelfilename)["beta"]
	lbfeatures = preprocesslearningbenchmarkfeatures()
	features, featureinfo, featurenums = [], [], []
else
	beta, features, featuresfor, featurenums, featureinfo = getmlfeatures(mlmodelfilename)
end

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
numpartitions, partitions, partitioninfo, globalpartitionid, partitionsolvetime = decomposeproblem(stationsperpartition, partitionobjective, beta, features, featureinfo, featurenums)
globalpartition = partitioninfo[globalpartitionid]
println("partitionsolvetime = $partitionsolvetime")

#Report 
writeglobalsolutionoutputs_init(globalsolutionfilename)
writeglobalsolutionoutputs_partitioning(partitionsolvetime)

#Enumerate windows for the global partition
windows_global, windowsduring_global, windowidlookup_global, windowscontaining_global = enumeratesubproblemwindows(globalpartition, 1, subproblemtimelength)
if shortmethodname != "learnbench"
	windowsynergy_global = preprocesswindowsynergies(windows_global, windowidlookup_global, featureinfo, features, beta, featurenums)
end

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

#Solve each partition
counter = 6
globalstarttime = time()
for s in 6:numpartitions

	println("===== PARTITION $s =====")

	#Get partition info
	currpartition = partitioninfo[s]

	#Find an initial solution
	currsol = partitionsolution[s] #createemptysolution(currpartition.orders, currpartition.pods, currpartition.podswith, currpartition.workstations)
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
	if shortmethodname != "learnbench"
		windowsynergy = preprocesswindowsynergies(windows, windowidlookup, featureinfo, features, beta, featurenums)
	end

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
		elseif shortmethodname == "learnbench"
			sp_winid, sp_orders, sp_window, sp_pods, sp_itemson, sp_items, tabulist, predicted_obj = selectlearningbenchmarksubproblem(currpartition, currsol, windows, windowidlookup, tabulist, sp_iter)
		end
		spselectiontime = time() - spselectionstarttime

		#Build subproblem structure
		sp = constructsubproblem(currpartition, sp_orders, sp_window, sp_pods, sp_itemson, sp_items, currsol)

		#Re-optimize subproblem
		sp_obj, sp_solvetime, h_sp, y_sp, z_sp, f_sp, g_sp, v_sp, feasibleflag_sp, sp_buildtime = reoptimizesubproblem(sp, currsol, currpartition, nocongestion_flag)
		println("Re-opt time = ", sp_solvetime, " seconds")

		#Update solution
		currsol = updatesolution(sp, currsol, currpartition, sp_obj, h_sp, y_sp, v_sp)
		updatecongestionfeatures()
		lastoptimizeddifference = updatelastoptimizeddifference(lastoptimizeddifference, tabulist, sp_winid, windows, predicted_obj, sp_obj)
		iterationtime = time() - iterationstarttime
		updatesolvemetrics(s, iterationtime, sp_solvetime, spselectiontime)

		#Write solution metrics
		writeglobalsolutionoutputs_iter(sp_iter, 0, iterationtime, spselectiontime, sp_solvetime, globalsolutionfilename, s, currpartition, currsol)

		#----------- REPORTING -----------#

		#Subproblem statistics
		if subproblemstatsreporting_flag == 1
			itemspicked = sum(sum(length(currsol.itempodpicklist[w,t]) for t in sp.times) for w in sp.workstations)
			podspicked = sum(sum(length(currsol.podsworkedat[w,t]) for t in sp.times) for w in sp.workstations)
			sp_util = (itemprocesstime*itemspicked + podprocesstime*podspicked) / sum(sum(tstep for t in sp.times) for w in sp.workstations) 
			sp_pileon = itemspicked/podspicked
			df = DataFrame(row_id=row_id, instance_id=instance_id, warehouse_id=warehouse_id, targetnumworkstations=targetnumworkstations, subproblemtimelength=subproblemtimelength, targetnumorders=targetnumorders, targetnumpods=targetnumpods, targetnumitems=targetnumitems, solvetime = sp_solvetime, buildtime = sp_buildtime, sp_obj = sp_obj, new_util = sp_util, new_pileon = sp_pileon, itemspicked=itemspicked, podspicked=podspicked)
			if sp_iter == 1
				CSV.write(subproblemstatsreportingfilename, df)
			else
				CSV.write(subproblemstatsreportingfilename, df, append=true)
			end
		end

		#Visualize new solution
		if visualizationflag == 1
			workstationviz(string(visualizationfolder, "/station_partition", s,"_iter", sp_iter,".png"), currpartition, currsol)
		end

		#For debugging: check consistency of current solution
		if debugmode == 1
			checksolution(currpartition, currsol, debugprintstatements)
		end

		println("Partition throughput = ", sum(length(currsol.itempodpicklist[w,t]) for w in currpartition.workstations, t in times))

		#Time limit termination (2 hours)
		#if time() - globalstarttime >= 60*60*2
		#	break
		#end

		global counter += 1

	end

	#Update global solution
	updateglobalsolution(currpartition, currsol)
	partitionsolution[s] = currsol

end
#-----------------------------------------------------------------------------------#

include("scripts/decomposition/writeglobalsolutionoutputs.jl")
writeglobalsolutionoutputs(globalsolutionfilename, solvemetrics)
# writepickdistrib(string(outputfolder,"/pickdistrib.csv"), globalsolution)
# writedistancedistrib(string(outputfolder,"/distdistrib.csv"), globalsolution)

#-----------------------------------------------------------------------------------#

if ordergraphreporting_flag == 1
	include("scripts/figures/ordergraph.jl")
	writeordergraph(ordergraphreportingfilename)
end

#-----------------------------------------------------------------------------------#
#=
include("scripts/routing/oldcongestionfunctions.jl")
include("scripts/routing/createroutingnetwork.jl")
include("scripts/routing/saveandloadassignments.jl")
savetaskassignments(string(outputfolder, "/assignments.jld2"))

include("run_routing.jl")
=#
#-----------------------------------------------------------------------------------#

#congestionanalysis(nocongestion_flag)
#include("scripts/visualizations/workstationviz.jl")
#workstationviz(string(visualizationfolder,"/station.png"), partitioninfo[1], partitionsolution[1])
#workstationviz_three(string(visualizationfolder,"/station_three.png"), partitioninfo[1], partitionsolution[1])

println("Done!")

#-----------------------------------------------------------------------------------#

#=
function findpath(n1)

	nodecosts = [999999 for n in 1:numnodes]
	nodecosts[n1] = 0
	for a in sort(1:numarcs, by=x->nodelookup[arclookup[x][1]][2])
		nstart, nend = arclookup[a]
		if nodecosts[nend] > nodecosts[nstart] - 1e-10
			nodecosts[nend] = nodecosts[nstart]
		end
	end

	return nodecosts

end

include("scripts/visualizations/timespacenetworkviz.jl")
numlocs = maximum(workstations)
for p in partitioninfo[1].pods[1:10]
	#wtowarcs = [a for a in 1:numarcs if (nodelookup[arclookup[a][1]][1] in workstations) & (nodelookup[arclookup[a][2]][1] in workstations) & (nodelookup[arclookup[a][1]][1] != nodelookup[arclookup[a][2]][1])]
	tsnarcs = [a for a in podarcset[p] if nodelookup[arclookup[a][2]][2] <= 360]
	pathexists = findpath(nodes[podstorageloc[p],0])
	for a in podarcset[p]
		if (arclookup[a][1] > numnodes) || (pathexists[arclookup[a][1]] == 999999)
			remove!(tsnarcs, a)
		end
	end
	usedarcs = union([arcs[nodes[1,t],nodes[1,t+tstep]] for t in 0:tstep:90]
	,[arcs[nodes[1,120],nodes[14,120+arclength[1,14]]], arcs[nodes[14,120+arclength[1,14]],nodes[11,120+arclength[1,14]+arclength[14,11]]]]
	,[arcs[nodes[11,t],nodes[11,t+tstep]] for t in 210:tstep:360-tstep])

	arclistlist = [tsnarcs, usedarcs]
	colorlist = [(150,150,150), (132,14,14)]
	thicknesslist = [5,10]
	dashlist = ["solid", "solid"]
	fractlist = [0,0]
	timespacenetwork(string(visualizationfolder, "/tsn_pod", p,".png"), arclistlist, colorlist, thicknesslist, dashlist, fractlist, 4000, 2000)
end

for p in partitioninfo[1].pods
	p=48
	wtowarcs = [a for a in 1:numarcs if (nodelookup[arclookup[a][1]][1] in workstations) & (nodelookup[arclookup[a][2]][1] in workstations) & (nodelookup[arclookup[a][1]][1] != nodelookup[arclookup[a][2]][1])]
	tsnarcs = setdiff(podarcset[p], wtowarcs)
	usedarcs = [a for a in podarcset[p] if partitionsolution[1].y[p,a] > 1e-4]

	arclistlist = [tsnarcs, wtowarcs, usedarcs]
	colorlist = [(180,180,180), (180,180,180), (132,14,14)]
	thicknesslist = [3,3,10]
	dashlist = ["solid", "solid", "solid"]
	fractlist = [0,0,0]
	timespacenetwork(string(visualizationfolder, "/tsn_pod", p,"_path.png"), arclistlist, colorlist, thicknesslist, dashlist, fractlist, 4000, 2000)
end

#Part 2
include("scripts/visualizations/timespacenetworkviz.jl")
p=1
numlocs = maximum(workstations)
tsnarcs = [a for a in podarcset[p] if nodelookup[arclookup[a][2]][2] <= 360]
pathexists = findpath(nodes[podstorageloc[p],0])
for a in podarcset[p]
	if (arclookup[a][1] > numnodes) || (pathexists[arclookup[a][1]] == 999999)
		remove!(tsnarcs, a)
	end
end
usedarcs = union([arcs[nodes[1,t],nodes[1,t+tstep]] for t in 0:tstep:90]
,[arcs[nodes[1,120],nodes[14,120+arclength[1,14]]], arcs[nodes[14,120+arclength[1,14]],nodes[11,120+arclength[1,14]+arclength[14,11]]]]
,[arcs[nodes[11,t],nodes[11,t+tstep]] for t in 210:tstep:360-tstep])

arclistlist = [tsnarcs, usedarcs]
colorlist = [(180,180,180), (132,14,14)]
thicknesslist = [5,10]
dashlist = ["solid", "solid"]
fractlist = [0,0]
timespacenetwork(string(visualizationfolder, "/tsn_anystorage_pod", p,".png"), arclistlist, colorlist, thicknesslist, dashlist, fractlist, 4000, 2000)

#Part 1
include("scripts/visualizations/timespacenetworkviz.jl")
p=1
numlocs = maximum(workstations)
tsnarcs = [a for a in podarcset[p] if (nodelookup[arclookup[a][2]][2] <= 360) & (nodelookup[arclookup[a][1]][1] in [1,13,14]) & (nodelookup[arclookup[a][2]][1] in [1,13,14])]
pathexists = findpath(nodes[podstorageloc[p],0])
for a in podarcset[p]
	if (arclookup[a][1] > numnodes) || (pathexists[arclookup[a][1]] == 999999)
		remove!(tsnarcs, a)
	end
end
usedarcs = union([arcs[nodes[1,t],nodes[1,t+tstep]] for t in 0:tstep:90]
,[arcs[nodes[1,120],nodes[14,120+arclength[1,14]]], arcs[nodes[14,120+arclength[1,14]],nodes[1,120+arclength[1,14]+arclength[14,1]]]]
,[arcs[nodes[1,t],nodes[1,t+tstep]] for t in 120+arclength[1,14]+arclength[14,1]:tstep:360-tstep])

arclistlist = [tsnarcs, usedarcs]
colorlist = [(180,180,180), (132,14,14)]
thicknesslist = [5,10]
dashlist = ["solid", "solid"]
fractlist = [0,0]
timespacenetwork(string(visualizationfolder, "/tsn_onlystorage_pod", p,".png"), arclistlist, colorlist, thicknesslist, dashlist, fractlist, 4000, 2000)


currsol = partitionsolution[1]
currpartition = partitioninfo[1]

totalpoddist, totalpodpicks = Dict(), Dict()
for p in currpartition.pods
	currtrippicks, currtripdist = 0, 0
	for a in currsol.ypath[p]
		l1,l2 = nodelookup[arclookup[a][1]][1], nodelookup[arclookup[a][2]][1]
		currtripdist += warehousedistance[l1,l2]
	end
	for w in workstations, t in times, (m,i,p2) in currsol.itempodpicklist[w,t]
		if p == p2
			currtrippicks += 1
		end
	end
	totalpoddist[p] = currtripdist
	totalpodpicks[p] = currtrippicks
end

df = (pod=[p for p in pods if totalpodpicks[p] >= 1], itempicks=[totalpodpicks[p] for p in pods if totalpodpicks[p] >= 1], disttraveled=[totalpoddist[p] for p in pods if totalpodpicks[p] >= 1])
CSV.write("figures/histograms/distdistrib_new.csv", df)


currsol = partitionsolution[1]
currpartition = partitioninfo[1]


totalpoddist, totalpodpicks = Dict(), Dict()
for p in currpartition.pods
	tripcounter = 1
	currtrippicks, currtripdist = 0, 0

	for a in currsol.ypath[p]
		l1,t1 = nodelookup[arclookup[a][1]]
		l2,t2 = nodelookup[arclookup[a][2]]
		currtripdist += warehousedistance[l1,l2]
		for w in workstations, t in max(0,t1):tstep:min(t2-tstep,horizon), (m,i,p2) in currsol.itempodpicklist[w,t]
			if p == p2
				currtrippicks += 1
			end
		end
		if (l2 == podstorageloc[p]) & !(l1 == podstorageloc[p])
			totalpoddist[p,tripcounter] = currtripdist
			totalpodpicks[p,tripcounter] = currtrippicks
			tripcounter += 1
			currtrippicks, currtripdist = 0, 0
		end
	end
	if currtrippicks > 1e-4
		totalpoddist[p,tripcounter] = currtripdist
		totalpodpicks[p,tripcounter] = currtrippicks
	end
end

allkeys = keys(totalpoddist)
df = (pod=[p for (p,ti) in allkeys if totalpodpicks[p,ti] >= 1], itempicks=[totalpodpicks[p,ti] for (p,ti) in allkeys if totalpodpicks[p,ti] >= 1], disttraveled=[totalpoddist[p,ti] for (p,ti) in allkeys if totalpodpicks[p,ti] >= 1])
CSV.write("figures/histograms/distdistrib_new.csv", df)

totalitempicks = Dict()
totalorderpicks = Dict()
for p in currpartition.pods
	pickeditems, pickedorders = [], []
	for w in workstations, t in times, (m,i,p2) in currsol.itempodpicklist[w,t]
		if p == p2
			push!(pickeditems, i)
			push!(pickedorders, m)
		end
	end
	totalitempicks[p] = length(pickeditems)
	totalorderpicks[p] = length(unique(pickedorders))
end


fullpodlist = [p for p in pods if totalitempicks[p] >= 1]
fulltypelist = ["Item picks" for p in pods if totalitempicks[p] >= 1]
fullpicklist = [totalitempicks[p] for p in pods if totalitempicks[p] >= 1]
for p in pods
	if totalitempicks[p] >= 1
		push!(fullpodlist, p)
		push!(fulltypelist, "Order picks")
		push!(fullpicklist, totalorderpicks[p])
	end
end


df = (pod=fullpodlist, picktype=fulltypelist, picks=fullpicklist)
CSV.write("figures/histograms/pickdistrib_new.csv", df)

=#

