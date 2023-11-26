
include("scripts/lsns/subproblemselection/findassignedorders.jl")

#-----------------------------------------------------------------------------------------------------#

function filloutrandomorders1(partition, sp_window, currsol, assignedorders, unassignedorders, numadditionalorders)

	sp_itemson = Dict()
	for m in assignedorders
		sp_itemson[m] = []
		for i in itemson[m]
			#if (sum(sum(currsol.h[m,i,p,w,sp_window.tend] for p in partition.podswith[i]) for w in sp_window.workstations) > 0.01) && (sp_window.tstart >= 0) && (sum(sum(currsol.h[m,i,p,w,sp_window.tstart-tstep] for p in partition.podswith[i]) for w in sp_window.workstations) < 0.01)
			#	push!(sp_itemson[m], i)
			if checkiteminpicklist(m, i, max(0,sp_window.tstart), min(horizon,sp_window.tend), sp_window.workstations, currsol)
				push!(sp_itemson[m], i)
				#println("FOUND AN ISSUE = $m, $i")
			elseif (sp_window.tstart == -30) && (sum(sum(currsol.h[m,i,p,w,sp_window.tend] for p in partition.podswith[i]) for w in sp_window.workstations) > 0.01)
				push!(sp_itemson[m], i)
			end
		end
	end

	additionalorders = unassignedorders[randperm(length(unassignedorders))][1:numadditionalorders]
	for m in additionalorders
		sp_itemson[m] = copy(itemson[m])
	end

	sp_orders = union(assignedorders, additionalorders)

	sp_items = []
	for m in sp_orders
		sp_items = union(sp_items, sp_itemson[m])
	end

	return sp_orders, sp_itemson, sp_items

end

#-----------------------------------------------------------------------------------------------------#

function findrandompods1(partition, targetnumpods, currsol, sp_orders, sp_window, sp_items)
	
	sp_pods = []
	for w in sp_window.workstations, t in sp_window.times
		sp_pods = union(sp_pods, currsol.podsworkedat[w,t])
	end

	potentialpods = setdiff(partition.pods, sp_pods) 
	numadditionalpods = max(targetnumpods - length(sp_pods), 0)
	extrapods = potentialpods[randperm(length(potentialpods))][1:numadditionalpods]
	sp_pods = union(sp_pods, extrapods)

	return sp_pods

end

#-----------------------------------------------------------------------------------------------------#

function findnecessarypods1(partition, targetnumpods, currsol, sp_orders, sp_window, sp_items)
	
	sp_pods = []
	for w in sp_window.workstations, t in sp_window.times
		sp_pods = union(sp_pods, currsol.podsworkedat[w,t])
	end

	usefulpods = []
	for m in sp_orders, i in itemson[m]
		usefulpods = union(usefulpods, partition.podswith[i])
	end

	potentialpods = setdiff(usefulpods, sp_pods) 
	numadditionalpods = max(targetnumpods - length(sp_pods), 0)
	extrapods = potentialpods[randperm(length(potentialpods))][1:numadditionalpods]
	sp_pods = union(sp_pods, extrapods)

	return sp_pods

end

#-----------------------------------------------------------------------------------------------------#

function selectrandomsubproblem(partition, windows, currsol, targetnumorders, targetnumpods)

	#Select the subproblem with the synergy model
	sp_window = rand(windows)

	#Get the orders
	assignedorders, unassignedorders, numadditionalorders = findassignedorders(partition, sp_window, currsol, targetnumorders)
	sp_orders, sp_itemson, sp_items = filloutrandomorders1(partition, sp_window, currsol, assignedorders, unassignedorders, numadditionalorders)
	
	#Get the pods
	sp_pods = findnecessarypods1(partition, targetnumpods, currsol, sp_orders, sp_window, sp_items)

	return sp_orders, sp_window, sp_pods, sp_itemson, sp_items

end
