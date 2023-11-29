
function findsubproblemnetwork(sp_workstations, sp_storagelocs, sp_tstart, sp_tend)

	tb_arcset, ap_arcset, sp_nodeset = [], [], []

	sp_times = max(0,sp_tstart):tstep:min(horizon,sp_tend)
	sp_alltimes = sp_tstart:tstep:sp_tend

	for t in sp_tstart:tstep:sp_tend, l in union(sp_workstations, sp_storagelocs)
		push!(sp_nodeset, extendednodes[l,t])
	end

	for n in sp_nodeset, a in A_plus[n]
		endloc, endtime = nodelookup[arclookup[a][2]]
		if (endtime <= sp_tend) & (endloc in union(sp_workstations, sp_storagelocs))
			push!(ap_arcset, a)
		end
	end

	for t in sp_tstart:tstep:sp_tend, l in union(workstations, storagelocs)
		n = extendednodes[l,t]
		for a in A_plus[n]
			if nodelookup[arclookup[a][2]][2] <= sp_tend
				push!(tb_arcset, a)
			end
		end
	end

	return tb_arcset, ap_arcset, sp_nodeset, sp_times, sp_alltimes

end

#-----------------------------------------------------------------------------------------------------#

function findsubproblempodsupply(partition, sp_times, sp_workstations, currsol, sp_arcset)

	sp_podsupply = Dict()

	for p in partition.pods
		for l in union(sp_workstations, podstorageloc[p]), t in sp_times
			sp_podsupply[p, extendednodes[l,t]] = 0
		end
	end

	for p in partition.pods
		for a in intersect(Set(currsol.ypath[p]), sp_arcset)
			snode, enode = arclookup[a]
			sp_podsupply[p, snode] -= 1
			sp_podsupply[p, enode] += 1
		end
	end

#	ERROR: LoadError: KeyError: key (1, 4341) not found

  	return sp_podsupply

end

#-----------------------------------------------------------------------------------------------------#

function findknownpodarcs(partition, currsol, sp_workstations, sp_tstart, sp_tend, sp_arcset)

	y_known, sp_podarcset, sp_podarcset_cong = Dict(), Dict(), Dict()

	for p in partition.pods
		sp_podarcset[p] = []
		sp_podarcset_cong[p] = []
		for a in intersect(podarcset[p], sp_arcset)
			push!(sp_podarcset[p], a)
		end
		for a in intersect(podarcset[p], sp_arcset, 1:numarcs)
			push!(sp_podarcset_cong[p], a)
		end
	end

	for p in partition.pods, w in sp_workstations, t in sp_tstart:tstep:sp_tend
		n = extendednodes[w,t]
		if setdiff(A_plus_p[p,n], union(A_queues, sp_arcset)) == []
			y_known[p,n] = 0
		else
			y_known[p,n] = sum(currsol.y[p,a] for a in setdiff(A_plus_p[p,n], union(A_queues, sp_arcset)))
		end
	end

	return y_known, sp_podarcset, sp_podarcset_cong

end

#-----------------------------------------------------------------------------------------------------#

function findsubproblemremaininginventory(partition, currsol, sp_workstations, sp_tstart, sp_tend)

	remaininginventory = deepcopy(inventory)

	for w in setdiff(partition.workstations, sp_workstations), t in times
		for item in currsol.itempodpicklist[w,t]
			m, i, p = item
			remaininginventory[i,p] -= 1
		end
	end

	for w in sp_workstations, t in setdiff(times, [t2 for t2 in sp_tstart:tstep:sp_tend])
		for item in currsol.itempodpicklist[w,t]
			m, i, p = item
			remaininginventory[i,p] -= 1
		end
	end

	return remaininginventory

end

#-----------------------------------------------------------------------------------------------------#

function getadditionalordersets(sp_orders, sp_itemson, sp_items)

	sp_orderswith = Dict()
	for i in sp_items
		sp_orderswith[i] = []
	end
	for m in sp_orders, i in sp_itemson[m]
		push!(sp_orderswith[i], m)
	end

	#Identify orders who are being worked both inside AND outside of this time window
	#We MUST complete these orders in the subproblem to ensure global feasibility
	sp_ordersinprogress = []
	for m in sp_orders
		if length(sp_itemson[m]) != length(itemson[m])
			push!(sp_ordersinprogress, m)
		end
	end

	#Use this version for julia 1.5 and lower:
	sp_numitems = 0
	for m in sp_orders
		sp_numitems += length(sp_itemson[m])
	end
	#Use this version for julia 1.6 and higher:		
	#sp_numitems = sum(length(sp_itemson[m]) for m in sp_orders; init=0)

	return sp_ordersinprogress, sp_orderswith, sp_numitems

end

#-----------------------------------------------------------------------------------------------------#

function findsubproblemambientcongestion(partition, currcong, sp_pods)

	#=
	ambientcongestion, emptycongestion = Dict(), Dict()

	excludedarcs = intersect(setdiff(tb_arcset, sp_arcset), 1:numarcs)
	for l in partition.intersections, t in max(0,sp_tstart):congestiontstep:min(horizon,sp_tend)
		ambientcongestion[l,t] = 0 
		emptycongestion[l,t] = 0 
		for p in partition.pods, a in intersect(partitionfloor.congestionarcs[l,t], podarcset[p], excludedarcs) 
			ambientcongestion[l,t] += partitionfloor.congestioncontribution[a,l,t] * currsol.y[p,a] 
		end
		#ambientcongestion[l,t] = sum(sum(partitionfloor.congestioncontribution[a,l,t] * currsol.y[p,a] for a in intersect(partitionfloor.congestionarcs[l,t], podarcset[p], excludedarcs)) for p in pods if intersect(partitionfloor.congestionarcs[l,t], podarcset[p], excludedarcs) != [])
	end
	=#

	ambientcongestion = sum(currcong[p] for p in setdiff(pods, sp_pods))
	emptycongestion = zeros(length(partition.intersections), length(0:congestiontstep:horizon))

	return ambientcongestion, emptycongestion

end

#-----------------------------------------------------------------------------------------------------#

function subprobleminventory(remaininginventory, podswith, sp_pods)
	
	sp_remaininginventory = Dict()
	for k in keys(remaininginventory)
		if k[2] in sp_pods
			sp_remaininginventory[k] = copy(remaininginventory[k])
		end
	end

	sp_podswith = deepcopy(podswith)
	for i in keys(podswith)
		counter = 1
		index = 1
		while counter <= length(podswith[i])
			if !(sp_podswith[i][index] in sp_pods)
				deleteat!(sp_podswith[i], index)
			else
				index = index+1
			end
			counter=counter+1
		end
	end

	return sp_remaininginventory, sp_podswith
end

#-----------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------#

function constructsubproblem(partition, sp_orders, sp_window, sp_pods, sp_itemson, sp_items, currsol)

	sp_workstations = sp_window.workstations
	sp_tstart = sp_window.tstart
	sp_tend = sp_window.tend

	st1 = time()
	tb_arcset, sp_arcset, sp_nodeset, sp_times, sp_alltimes = findsubproblemnetwork(sp_workstations, partition.storagelocs, sp_tstart, sp_tend)
	println("  Time 1 = ", time()-st1)
	st2 = time()
	sp_podsupply = findsubproblempodsupply(partition, sp_alltimes, sp_workstations, currsol, sp_arcset)
	println("  Time 2 = ", time()-st2)
	st3 = time()
	y_known, sp_podarcset, sp_podarcset_cong = findknownpodarcs(partition, currsol, sp_workstations, sp_tstart, sp_tend, sp_arcset)
	println("  Time 3 = ", time()-st3)
	st4 = time()
	remaininginventory = findsubproblemremaininginventory(partition, currsol, sp_workstations, sp_tstart, sp_tend)
	println("  Time 4 = ", time()-st4)
	st5 = time()
	sp_remaininginventory, sp_podswith = subprobleminventory(remaininginventory, podswith, sp_pods)
	println("  Time 5 = ", time()-st5)
	st6 = time()
	sp_ordersinprogress, sp_orderswith, sp_numitems = getadditionalordersets(sp_orders, sp_itemson, sp_items)
	println("  Time 6 = ", time()-st6)
	st7 = time()
	ambientcongestion, emptycongestion = findsubproblemambientcongestion(partition, currcong, sp_pods)
	println("  Time 7 = ", time()-st7)

	sp = (orders=sp_orders, workstations=sp_workstations, pods=sp_pods, tstart=sp_tstart, tend=sp_tend, partition=partition,
		timeblockarcset=tb_arcset, arcset=sp_arcset, nodeset=sp_nodeset, times=sp_times, alltimes=sp_alltimes, podsupply=sp_podsupply,
		y_known=y_known, podarcset=sp_podarcset, sp_podarcset_cong=sp_podarcset_cong, remaininginventory=remaininginventory, podswith=sp_podswith, ambientcongestion=ambientcongestion,
		itemson=sp_itemson, orderswith=sp_orderswith, items=sp_items, ordersinprogress=sp_ordersinprogress
		)

	return sp

end
