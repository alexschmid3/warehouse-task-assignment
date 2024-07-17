
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
include("scripts/lsns/createemptysolution.jl")
include("scripts/lsns/reoptimizesubproblem.jl")
include("scripts/lsns/updatesolution.jl")
include("scripts/lsns/modelfeatures/getlearningbenchmarkfeatures.jl")
include("scripts/lsns/modelfeatures/preprocesslearningbenchmarkfeatures.jl")
include("scripts/lsns/modelfeatures/getmlfeatures.jl")
include("scripts/lsns/subproblemselection/selectrandomsubproblem.jl")
include("scripts/lsns/subproblemselection/selectsynergisticsubproblem.jl")
include("scripts/lsns/subproblemselection/selectlearnthenoptimizesubproblem.jl")
include("scripts/lsns/subproblemselection/selectlearningbenchmarksubproblem.jl")
include("scripts/test/checksolution.jl")
include("scripts/figures/histograms.jl")
include("scripts/training/datageneration/findimpactordersandpods.jl")
include("scripts/training/datageneration/getfeatures.jl")
include("scripts/training/datageneration/staticfeatures.jl")
include("scripts/training/datageneration/savesubproblembuild.jl")

println("Scripts imported")

#-----------------------------------------------------------------------------------#

#Debugging mode
visualizationflag = 0			# 1 --> Produce visualizations of warehouse and workstation solutions
debugprintstatements = 0 		# When a bug is found, set this to 1 for next run to get more detailed error info (warning: it's a lot of print statements, one for each unit test)
debugmode = 0					# 1 --> will perform solution consistency unit tests at each LSNS iteration, use for debugging when changing the code, 0 --> no solution checks, use for computational experiments

#Initialize Gurobi
const GRB_ENV = Gurobi.Env()

# Select the instancecd
row_id = ifelse(length(ARGS) > 0, parse(Int, ARGS[1]), 1) # (for cluster submissions)
warehouseparamsfilename = "data/warehouse_sizes_and_capacities.csv"
instanceparamsfilename = "data/train_instance_parameters.csv"
methodparamsfilename = "data/train_run_parameters.csv"
warehouseparms = CSV.read(warehouseparamsfilename, DataFrame)
instanceparms = CSV.read(instanceparamsfilename, DataFrame)
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
stationtostation_flag = methodparms[row_id, 11]
targetnumworkstations = methodparms[row_id, 12]
mlmodelfilename = ""
minnumpods = targetnumpods - 10
minnumorders = targetnumorders - 10
minnumitems = targetnumitems - 10
numlocalintersections = 6
numhyperlocalintersections = 3 
maxworkstationspersubproblem = 2
shortmethodname, subproblembudget = parsemethodname(methodname)
anystoragelocation_flag = 0

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
partitionobjective = "none"

#Training data
dynamicmlpass = 0
staticmlpass = 0
traindatafolder = string("trainingdata/mainmodel_wh",warehouse_id)
dynamicinstanceoutputfilename = string(traindatafolder,"/dynamic/features_wh", warehouse_id, "_pass", dynamicmlpass, "_instance", instance_id, "_run", run_id,".jld2")
staticinstanceoutputfilename = string(traindatafolder,"/static/features_wh", warehouse_id, "_pass", staticmlpass, "_instance", instance_id, "_run", run_id, ".csv")

println("Parameters read")

#Files
outputfolder = string("outputs/mainmodeltrainingrun", run_id,"_", today())
globalsolutionfilename = string(outputfolder, "/output.csv")
if !(isdir(outputfolder))
	mkdir(outputfolder)
end
if !(isdir(traindatafolder))
	mkdir(traindatafolder)
	mkdir(string(traindatafolder,"/dynamic"))
	mkdir(string(traindatafolder,"/static"))
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

#ML features
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
#windows_global, windowsduring_global, windowidlookup_global, windowscontaining_global = enumeratesubproblemwindows(globalpartition, 1, subproblemtimelength)
#if shortmethodname != "learnbench"
#	windowsynergy_global = preprocesswindowsynergies(windows_global, windowidlookup_global, featureinfo, features, beta, featurenums)
#end

#-----------------------------------------------------------------------------------#

#Static
subproblemfeatures_df = createemptydataframe()

#Dynamic
subproblemarray = []
x_k, y_k, z_k, v_k, q_k, obj_k = Dict(), Dict(), Dict(), Dict(), Dict(), Dict()
k_star, best_improvement = 0, 0

#Algorithm congtol 
windowrecentlyreoptimized = zeros(numsubproblemsevaluated)
windowreoptimizedcount = zeros(numsubproblemsevaluated)
mostrecentwindowimprovement = [1 for j in 1:numsubproblemsevaluated]
recentlyoptimizedwindows = []

#Static
orderordercompat = calculateordercompatibility()
orderpodcompat = calculateorderpodcompatibility()
centrality, horizonrelativetime, avgdisttostations, ordersize, oneitemflag, iteminventory, itemnumpods = preprocessentityfeatures_static()
itemoverlap, itemoverlappct, overlapitemavginv, overlapitemavgaltpods, pwdistance, xsteps, ysteps = preprocesssynergyfeatures_static()

#Dynamic
allpartition = (workstations=workstations, orders=orders, pods=pods)
#instance_features = getinstancefeatures(allpartition, currsol)
w_features, t_features, p_features, m_features, i_features = getentityfeatures()
mp_features, pw_features, pwt_features, wt_features = getsynergyfeatures()
#currsol_features, openorders, openpods = getcurrsolfeatures()

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

#Solve each partition
counter = 1
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
	#if shortmethodname != "learnbench"
	#	windowsynergy = preprocesswindowsynergies(windows, windowidlookup, featureinfo, features, beta, featurenums)
	#end

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
			sp_winid, sp_orders, sp_window, sp_pods, sp_itemson, sp_items, assignedorders = selectrandomsubproblem(currpartition, windows, windowidlookup, currsol, targetnumorders, targetnumpods)
			predicted_obj = 0
		elseif methodname == "synergy"
			sp_winid, sp_orders, sp_window, sp_pods, sp_itemson, sp_items, tabulist = selectsynergisticsubproblem(currpartition, windows, currsol, targetnumorders, targetnumpods, rand(1:maxworkstationspersubproblem), tabulist)
			predicted_obj = 0
		elseif shortmethodname == "learnbench"
			sp_winid, sp_orders, sp_window, sp_pods, sp_itemson, sp_items, tabulist, predicted_obj = selectlearningbenchmarksubproblem(currpartition, currsol, windows, windowidlookup, tabulist, sp_iter)
		end
		spselectiontime = time() - spselectionstarttime
		println("Selection time = ", spselectiontime, " seconds")

		#Build subproblem structure
		sp = constructsubproblem(currpartition, sp_orders, sp_window, sp_pods, sp_itemson, sp_items, currsol)

		#Re-optimize subproblem
		sp_obj, sp_solvetime, h_sp, y_sp, z_sp, f_sp, g_sp, v_sp, feasibleflag_sp = reoptimizesubproblem(sp, currsol, currpartition, 0)
		println("Re-opt time = ", sp_solvetime, " seconds")

		#------------------ Add subproblem data ------------------#

		#Check objective improvement
		old_obj = sum(sum(length(currsol.itempodpicklist[w,t]) for w in sp.workstations) for t in max(0,sp.tstart):tstep:min(horizon,sp.tend))
		if (old_obj > 0) & (sp_obj > 0)
			println("Improved from ", old_obj, " to ", sp_obj, " = ", round(100*(sp_obj - old_obj) / old_obj, digits=2), "%")
		else
			println("Improved from ", old_obj, " to ", sp_obj, " = 0.00 %")
		end

		#Get the impacted orders and podsp
		if sp_obj > 0
			impact_orders, impact_pods = findimpactordersandpods(sp.workstations, sp.orders, sp.pods, sp.itemson, sp.times, sp.podswith, currsol.h, h_sp)
		else
			impact_orders, impact_pods = [], []
		end

        #Solve problem without congestion to detect congestion impacts
        nocong_obj, nocong_solvetime, h_nocong, y_nocong, z_nocong, f_nocong, g_nocong, v_nocong, feasibleflag_nocong = reoptimizesubproblem(sp, currsol, currpartition, 1)
        #Get the impacted orders and podsp
        if nocong_obj > 0
            nocong_impact_orders, nocong_impact_pods = findimpactordersandpods(sp.workstations, sp.orders, sp.pods, sp.itemson, sp.times, sp.podswith, currsol.h, h_nocong)
        else
            nocong_impact_orders, nocong_impact_pods = [], []
        end

        #Get the predictions and add to dataframe (STATIC)
        staticfeatures = getsubproblemfeatures(sp_iter, sp, 0, sp.workstations, sp.times, sp.tstart, sp.tend, sp.orders, sp.nodeset, sp.arcset, sp.itemson, sp.ambientcongestion, sp.y_known, assignedorders, sp.pods, currsol)
        outcomes = formatoutcomes(old_obj, sp_obj, sp_solvetime)
        newrow = hcat(staticfeatures, outcomes)
        push!(subproblemfeatures_df, newrow)

        #Record improvement (DYNAMIC)
        subprobleminfo = Dict("sp_opt" => 0, "sp_workstations" => sp.workstations, "sp_tstart" => sp.tstart, "sp_tend" => sp.tend, "sp_orders" => sp.orders, "sp_pods" => sp.pods, "sp_numitems" => length(sp.items), "impact_orders" => impact_orders, "impact_pods" => impact_pods, "old_obj" => old_obj, "new_obj" => sp_obj, "nocongestion_obj" => nocong_obj, "nocongestion_impact_orders" => nocong_impact_orders, "nocongestion_impact_pods" => nocong_impact_pods, "sp_solvetime" => sp_solvetime)
        push!(subproblemarray, subprobleminfo)
        savesubproblembuild(sp_iter, sp.workstations, sp.times, impact_orders, impact_pods, sp_obj, old_obj)

        #---------------------------------------------------------#
        
		#Update solution
		#spupdatestarttime = time()
		#currsol = updatesolution(sp, currsol, currpartition, sp_obj, h_sp, y_sp, v_sp)
		#lastoptimizeddifference = updatelastoptimizeddifference(lastoptimizeddifference, tabulist, sp_winid, windows, predicted_obj, sp_obj)
		#iterationtime = time() - iterationstarttime
		#updatesolvemetrics(s, iterationtime, sp_solvetime, spselectiontime)

		#Write solution metrics
		#writeglobalsolutionoutputs_iter(sp_iter, 0, iterationtime, spselectiontime, sp_solvetime, globalsolutionfilename, s, currpartition, currsol)

		#Visualize new solution
		#if visualizationflag == 1
		#	workstationviz(string(visualizationfolder, "/station_partition", s,"_iter", sp_iter,".png"), currpartition, currsol)
		#end

		#For debugging: check consistency of current solution
		if debugmode == 1
			checksolution(currpartition, currsol, debugprintstatements)
		end

		#println("Partition throughput = ", sum(length(currsol.itempodpicklist[w,t]) for w in currpartition.workstations, t in times))

		global counter += 1

	end

    #Save static
    CSV.write(staticinstanceoutputfilename, subproblemfeatures_df)

    #Save dynamic
    instance_features = getinstancefeatures(currpartition, currsol)
    save(dynamicinstanceoutputfilename, Dict("instance_features" => instance_features, "w_features" => w_features, "t_features" => t_features, "m_features" => m_features, "p_features" => p_features, "i_features" => i_features, "wt_features" => wt_features, "mp_features" => mp_features, "pw_features" => pw_features, "pwt_features" => pwt_features, "subproblemdata" => subproblemarray, "maps" => maps))

	#Update global solution
	#updateglobalsolution(currpartition, currsol)
	#partitionsolution[s] = currsol

end

#-----------------------------------------------------------------------------------#

println("Done!")