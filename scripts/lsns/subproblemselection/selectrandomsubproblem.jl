
function findassignedorders1(partition, sp_window, currsol, targetnumorders)

	assignedorders = []
	for w in sp_window.workstations, t in sp_window.times
		#Include all open orders, regardless of whether there are items picked in the window
		assignedorders = union(assignedorders, currsol.ordersopen[w,t])

		#Include all orders that have an item picked in the window, regardless of whether they remained open for a full period of the window
		#Note this is necessary to ensure single item orders are included
		for (m,i,p) in currsol.itempodpicklist[w,t]
			assignedorders = union(assignedorders, m)
		end
	end

	#Calculate how many more orders are needed to reach the target
	unassignedorders = [m for m in partition.orders if sum(sum(sum(currsol.h[m,i,p,w,horizon] for w in partition.workstations) for p in partition.podswith[i]) for i in itemson[m]) < 1e-4]
	unassignedorders = setdiff(unassignedorders, assignedorders)
	numadditionalorders = min(length(unassignedorders), max(0, targetnumorders - length(assignedorders)))
	
	return assignedorders, unassignedorders, numadditionalorders

end

#-----------------------------------------------------------------------------------------------------#

function checkiteminpicklist(m, i, tstart, tend, workstations, currsol)

	for w in workstations, t in tstart:tstep:tend, (m1,i1,p1) in currsol.itempodpicklist[w,t]
		if (m1 == m) & (i1 == i)
			return true
		end
	end

	return false

end

#-----------------------------------------------------------------------------------------------------#

function filloutrandomorders1(partition, sp_window, currsol, assignedorders, unassignedorders, numadditionalorders)

	sp_itemson = Dict()
	for m in assignedorders
		sp_itemson[m] = []
		for i in itemson[m]
			if (sum(sum(currsol.h[m,i,p,w,sp_window.tend] for p in partition.podswith[i]) for w in sp_window.workstations) > 0.01) && (sp_window.tstart >= 0) && (sum(sum(currsol.h[m,i,p,w,sp_window.tstart-tstep] for p in partition.podswith[i]) for w in sp_window.workstations) < 0.01)
				push!(sp_itemson[m], i)
			elseif checkiteminpicklist(m, i, max(0,sp_window.tstart), min(horizon,sp_window.tend), sp_window.workstations, currsol)
				push!(sp_itemson[m], i)
				println("FOUND AN ISSUE = $m, $i")
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

function selectrandomsubproblem(partition, windows, windowidlookup, currsol, targetnumorders, targetnumpods)

	#Select the subproblem with the synergy model
	sp_window = rand(windows)
	sp_winid = windowidlookup[sp_window]

	#Get the orders
	assignedorders, unassignedorders, numadditionalorders = findassignedorders1(partition, sp_window, currsol, targetnumorders)
	sp_orders, sp_itemson, sp_items = filloutrandomorders1(partition, sp_window, currsol, assignedorders, unassignedorders, numadditionalorders)
	
	#Get the pods
	sp_pods = findnecessarypods1(partition, targetnumpods, currsol, sp_orders, sp_window, sp_items)

	return sp_winid, sp_orders, sp_window, sp_pods, sp_itemson, sp_items

end
