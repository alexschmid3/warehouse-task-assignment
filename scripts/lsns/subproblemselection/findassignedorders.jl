
function findassignedorders(currpartition, sp_window, currsol, targetnumorders)

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
	unassignedorders = [m for m in currpartition.orders if sum(sum(sum(currsol.h[m,i,p,w,horizon] for w in currpartition.workstations) for p in currpartition.podswith[i]) for i in itemson[m]) < 1e-4]
	unassignedorders = setdiff(unassignedorders, assignedorders)
	numadditionalorders = min(length(unassignedorders), max(0, targetnumorders - length(assignedorders)))
	
	return assignedorders, unassignedorders, numadditionalorders

end

#-----------------------------------------------------------------------------------------------------#

function findallunassignedorders(currpartition, currsol)

	assignedorders = []
	for w in currpartition.workstations, t in times
		#Include all open orders, regardless of whether there are items picked in the window
		assignedorders = union(assignedorders, currsol.ordersopen[w,t])

		#Include all orders that have an item picked in the window, regardless of whether they remained open for a full period of the window
		#Note this is necessary to ensure single item orders are included
		for (m,i,p) in currsol.itempodpicklist[w,t]
			assignedorders = union(assignedorders, m)
		end
	end

	#Calculate how many more orders are needed to reach the target
	unassignedorders = [m for m in currpartition.orders if sum(sum(sum(currsol.h[m,i,p,w,horizon] for w in currpartition.workstations) for p in currpartition.podswith[i]) for i in itemson[m]) < 1e-4]
	unassignedorders = setdiff(unassignedorders, assignedorders)
	
	return assignedorders, unassignedorders

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

