
using NearestNeighbors 

#-----------------------------------HELPER FUNCTIONS------------------------------------# 

#Removing item from list
function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end

#-----------------------------------------------------------------------------------#

function arcdesc(a)

	println(nodelookup[arclookup[a][1]], " ==> ", nodelookup[arclookup[a][2]])

end

#-----------------------------------------------------------------------------------#

function getlocations(podstart, workloc)

	numstoragelocs, numworkstations = size(DataFrames.combine(DataFrames.groupby(podstart, [:start_x, :start_y]), nrow => :count))[1], size(workloc)[1]
	storagelocs = [l for l in 1:numstoragelocs]
	workstations = [l for l in numstoragelocs+1:numstoragelocs+numworkstations]
	
	xcoords = vcat(DataFrames.combine(DataFrames.groupby(podstart, [:start_x, :start_y]), nrow => :count)[!, 1], workloc[!, 2])
	ycoords = vcat(DataFrames.combine(DataFrames.groupby(podstart, [:start_x, :start_y]), nrow => :count)[!, 2], workloc[!, 3])
	loccoords = hcat(xcoords, ycoords)

	loclookup = Dict()
	for l in 1:size(loccoords)[1]
		loclookup[loccoords[l,1],loccoords[l,2]] = l
	end

	return storagelocs, workstations, loccoords, loclookup, numstoragelocs, numworkstations

end

#-----------------------------------------------------------------------------------#

function getintersections(horstreets, vertstreets, numstoragelocs, numworkstations)

	intersections, queueintersections = [], []
	intcoords, intlookup = Dict(), Dict()
	xcoords, ycoords = [], []

	#Storage space intersections
	index = numstoragelocs + numworkstations + 1
	for hs in 1:size(horstreets)[1], vs in 1:size(vertstreets)[1]
		inty = (horstreets[hs,1] + horstreets[hs,2]) / 2
		intx = (vertstreets[vs,1] + vertstreets[vs,2]) / 2
		push!(intersections, index)
		intcoords[index] = (intx, inty)
		intlookup[intx, inty] = index
		push!(xcoords, intx)
		push!(ycoords, inty)
		index += 1
	end

	#Workstation queue intersections
	for w in workstations
		wscoord = loccoords[w,:]
		inty = wscoord[2] + 2 
		intx = wscoord[1] + 7
		push!(intersections, index)
		push!(queueintersections, index)
		intcoords[index] = (intx, inty)
		intlookup[intx, inty] = index
		push!(xcoords, intx)
		push!(ycoords, inty)
		index += 1
	end

	intcoords_nn = hcat(xcoords, ycoords)

	return intersections, queueintersections, intcoords, intlookup, intcoords_nn
		
end

#-----------------------------------------------------------------------------------#

function maplocationaccesspoints(storagelocs, workstations, intcoords_nn, loccoords, intlookup)

	maploctointersection = Dict()

	kdtree = KDTree(transpose(intcoords_nn))

	for l in storagelocs
		idx, dists = knn(kdtree, loccoords[l,:], 1, true)
		maploctointersection[l] = intlookup[intcoords_nn[idx,:][1], intcoords_nn[idx,:][2]]
	end

	N = size(intcoords_nn)[1]
	kdtree_w = KDTree(transpose(intcoords_nn[N-numworkstations+1:N,:]))

	for w in workstations
		idx, dists = knn(kdtree_w, loccoords[w,:], 1, true)
		maploctointersection[w] = intlookup[intcoords_nn[idx[1]+N-numworkstations,:][1], intcoords_nn[idx[1]+N-numworkstations,:][2]]
	end

	return maploctointersection

end

#-----------------------------------------------------------------------------------#

function createtimespacenetwork(storagelocs, workstations, intersections)

	nodes, nodelookup, arcs, arclookup, A_plus, A_minus = Dict(), Dict(), Dict(), Dict(), Dict(), Dict()
	N_end, times = [], []

	index = 1
	for t in 0:tstep:horizon
		for l in storagelocs
			nodes[l,t] = index
			nodelookup[index] = (l,t)
			if t == horizon
				push!(N_end, index)
			end
			index += 1
		end	
		for w in workstations
			nodes[w,t] = index
			nodelookup[index] = (w,t)
			if t == horizon
				push!(N_end, index)
			end
			index += 1
		end	
		for l in intersections
			nodes[l,t] = index
			nodelookup[index] = (l,t)
			if t == horizon
				push!(N_end, index)
			end
			index += 1
		end	
	end
	
	numnodes = length(nodes)
	times = [t for t in 0:tstep:horizon]

	return numnodes, nodes, nodelookup, N_end, times

end

#-----------------------------------------------------------------------------------#

function extendtimespacenetwork(nodelookup)

	dummystarttime, dummyendtime = -tstep, horizon + tstep
	N_ext_before, N_ext_after = [], []
	extendednodes = deepcopy(nodes)
	index = numnodes + 1
	for t in [dummystarttime, dummyendtime], l in union(storagelocs, workstations)
		extendednodes[l,t] = index
		nodelookup[index] = (l,t)
		if t == dummystarttime
			push!(N_ext_before, index)
		else
			push!(N_ext_after, index)	
		end
		index += 1
	end

	extendedtimes = union(dummystarttime, times, dummyendtime)
	extendednumnodes = length(extendednodes)

	return extendednumnodes, extendednodes, nodelookup, N_ext_before, N_ext_after, extendedtimes, dummystarttime, dummyendtime

end

#-----------------------------------------------------------------------------------#

function getphysicalarcs(podspeed, loccoords)
	
	prearcs = []
	arclength = Dict()

	#Added some arc timing adjustments for consistency across various time discretizations
	if warehouseparamsfilename == "data/extensions/timedisc/warehouse_sizes_and_capacities.csv"
		for s in storagelocs, w in workstations
			wint = maploctointersection[w]
			manhattandist = abs(loccoords[s,1] - intcoords[wint][1]) +  abs(loccoords[s,2] - intcoords[wint][2])
			arclengthraw = manhattandist / podspeed
			leng = tstep * ceil(arclengthraw / tstep)
			push!(prearcs, (s, w, manhattandist, arclengthraw, leng))
			push!(prearcs, (w, s, manhattandist, arclengthraw, leng))
			arclength[s, w] = leng
			arclength[w, s] = leng
		end

		if stationtostation_flag == 1
			for w1 in workstations, w2 in [w for w in workstations if w < w1]
				manhattandist = 6 + abs(loccoords[w1,1] - loccoords[w2,1]) +  abs(loccoords[w1,2] - loccoords[w2,2])
				arclengthraw = manhattandist / podspeed
				leng = tstep * ceil(arclengthraw / tstep)
				push!(prearcs, (w1, w2, manhattandist, arclengthraw, leng))
				push!(prearcs, (w2, w1, manhattandist, arclengthraw, leng))
				arclength[w1, w2] = leng
				arclength[w2, w1] = leng
			end
		end

	#Original version of timing used for all other experiments
	else
		for s in storagelocs, w in workstations
			manhattandist = abs(loccoords[s,1] - loccoords[w,1]) +  abs(loccoords[s,2] - loccoords[w,2])
			arclengthraw = manhattandist / podspeed
			leng = tstep * ceil(arclengthraw / tstep)
			push!(prearcs, (s, w, manhattandist, arclengthraw, leng))
			push!(prearcs, (w, s, manhattandist, arclengthraw, leng))
			arclength[s, w] = leng
			arclength[w, s] = leng
		end

		if stationtostation_flag == 1
			for w1 in workstations, w2 in [w for w in workstations if w < w1]
				manhattandist = abs(loccoords[w1,1] - loccoords[w2,1]) +  abs(loccoords[w1,2] - loccoords[w2,2])
				arclengthraw = manhattandist / podspeed
				leng = tstep * ceil(arclengthraw / tstep)
				push!(prearcs, (w1, w2, manhattandist, arclengthraw, leng))
				push!(prearcs, (w2, w1, manhattandist, arclengthraw, leng))
				arclength[w1, w2] = leng
				arclength[w2, w1] = leng
			end
		end
	end

	return prearcs, arclength

end

#-----------------------------------------------------------------------------------#

function createtimespacearcs(prearcs, numnodes)

	arcs, arclookup, A_plus, A_minus = Dict(), Dict(), Dict(), Dict()
	A_space = []
	
	for node in 1:numnodes
		A_plus[node] = []
		A_minus[node] = []
	end

	prearcs_aug = deepcopy(prearcs)
	for l in union(storagelocs, workstations)
		push!(prearcs_aug, (l, l, 0, tstep, tstep))
	end

	index = 1
	for arc in prearcs_aug, t in 0:tstep:horizon-arc[5]
		startnode = nodes[arc[1],t]
		endnode = nodes[arc[2],t+arc[5]]
	
		arcs[(startnode,endnode)] = index
		arclookup[index] = (startnode,endnode)
		
		push!(A_plus[startnode], index)
		push!(A_minus[endnode], index)
		if arc in prearcs
			push!(A_space, index)
		end

		index += 1
	end

	numarcs = length(arcs)

	A_queues = []
	for w in workstations, t in setdiff(times,horizon) 
		A_queues = union(A_queues, [a for a in setdiff(A_plus[nodes[w,t]], A_space)])
	end

	return numarcs, arcs, arclookup, A_plus, A_minus, A_space, A_queues

end

#-----------------------------------------------------------------------------------#

function extendtimespacearcs(arclookup, A_space, A_plus, A_minus)

	extendedarcs = deepcopy(arcs)

	for n in union(N_ext_after, N_ext_before)
		A_plus[n] = []
		A_minus[n] = []
	end

	index = length(arcs) + 1
	for pa in prearcs, t2 in [t for t in 0:tstep:pa[5]]
		l1, l2 = pa[1], pa[2]
		t1 = dummystarttime
		arclookup[index] = (extendednodes[l1,t1], extendednodes[l2,t2])
		extendedarcs[extendednodes[l1,t1], extendednodes[l2,t2]] = index
		push!(A_plus[extendednodes[l1,t1]], index)
		push!(A_minus[extendednodes[l2,t2]], index)
		if l1 != l2
			push!(A_space, index)
		end
		index += 1
	end
	for pa in prearcs, t1 in [t for t in dummyendtime-pa[5]:tstep:horizon]
		l1, l2 = pa[1], pa[2]
		t2 = dummyendtime
		arclookup[index] = (extendednodes[l1,t1], extendednodes[l2,t2])
		extendedarcs[extendednodes[l1,t1], extendednodes[l2,t2]] = index
		push!(A_plus[extendednodes[l1,t1]], index)
		push!(A_minus[extendednodes[l2,t2]], index)
		if l1 != l2
			push!(A_space, index)
		end
		index += 1
	end
	for t1 in [dummystarttime, horizon], l in union(storagelocs, workstations)
		t2 = t1 + tstep
		arclookup[index] = (extendednodes[l,t1], extendednodes[l,t2])
		extendedarcs[extendednodes[l,t1], extendednodes[l,t2]] = index
		push!(A_plus[extendednodes[l,t1]], index)
		push!(A_minus[extendednodes[l,t2]], index)
		index += 1
	end

	extendednumarcs = length(extendedarcs)
	
	return extendednumarcs, extendedarcs, arclookup, A_plus, A_minus, A_space

end

#-----------------------------------------------------------------------------------#

function getpods(podstart)

	numpods = size(podstart)[1]
	allpods = [p for p in 1:numpods]
	podstartnode = [extendednodes[loclookup[podstart[p,2], podstart[p,3]], dummystarttime] for p in allpods]
	podstorageloc = [loclookup[podstart[p,2], podstart[p,3]] for p in allpods]

	return numpods, allpods, podstartnode, podstorageloc

end

#-----------------------------------------------------------------------------------#

function getinventory(inventoryinfo)

	allitems = [i for i in 1:num_unique_items]

	podswith = Dict()
	inventory = Dict()
	podstartinventory = Dict()
	for i in allitems
		podswith[i] = []
	end
	for p in allpods
		podstartinventory[p] = []
	end
	for row in 1:size(inventoryinfo)[1]
		push!(podswith[convert(Int64, inventoryinfo[row,2])], convert(Int64, inventoryinfo[row,1]))
		push!(podstartinventory[convert(Int64, inventoryinfo[row,1])], convert(Int64, inventoryinfo[row,2]))
		inventory[convert(Int64, inventoryinfo[row,2]), convert(Int64, inventoryinfo[row,1])] = inventoryinfo[row,3]
	end

	return podswith, allitems, inventory, podstartinventory	

end

#-----------------------------------------------------------------------------------#

function podarcsets(items, arclength)

	pods = copy(allpods)
	#pods = []
	#for i in items
	#	pods = union(pods, podswith[i])
	#end
	#pods = sort(pods)

	podnodeset = Dict()
	for p in pods
		podnodeset[p] = []
		if anystoragelocation_flag == 1
			for l in union(workstations, storagelocs), t in extendedtimes
				push!(podnodeset[p] , extendednodes[l,t])
			end
		else
			for l in union(workstations, podstorageloc[p]), t in extendedtimes
				push!(podnodeset[p] , extendednodes[l,t])
			end
		end
	end

	podarcset = Dict()
	arcpodset = Dict()
	A_minus_p, A_plus_p = Dict(), Dict()

	for p in pods
		podarcset[p] = []
	end
	for p in pods, n in 1:extendednumnodes
		A_plus_p[p,n] = []
		A_minus_p[p,n] = []
	end
	for a in 1:extendednumarcs
		arcpodset[a] = []
	end
	
	#Movement arcs
	if anystoragelocation_flag == 1
		for p in pods, w in workstations, s in storagelocs, t in 0:tstep:horizon-arclength[s, w]
			a_leave = arcs[nodes[s, t], nodes[w, t + arclength[s, w]]]
			a_return = arcs[nodes[w, t], nodes[s, t + arclength[w, s]]]
			
			push!(podarcset[p], a_leave)
			push!(podarcset[p], a_return)
			push!(arcpodset[a_leave], p)
			push!(arcpodset[a_return], p)
			push!(A_plus_p[p, nodes[s, t]], a_leave)
			push!(A_plus_p[p, nodes[w, t]], a_return)
			push!(A_minus_p[p, nodes[w, t + arclength[s, w]]], a_leave)
			push!(A_minus_p[p,nodes[s, t + arclength[w, s]]], a_return)
		end
	else
		for p in pods, w in workstations, t in 0:tstep:horizon-arclength[podstorageloc[p], w]
			s = podstorageloc[p]
			a_leave = arcs[nodes[s, t], nodes[w, t + arclength[s, w]]]
			a_return = arcs[nodes[w, t], nodes[s, t + arclength[w, s]]]
			
			push!(podarcset[p], a_leave)
			push!(podarcset[p], a_return)
			push!(arcpodset[a_leave], p)
			push!(arcpodset[a_return], p)
			push!(A_plus_p[p, nodes[s, t]], a_leave)
			push!(A_plus_p[p, nodes[w, t]], a_return)
			push!(A_minus_p[p, nodes[w, t + arclength[s, w]]], a_leave)
			push!(A_minus_p[p,nodes[s, t + arclength[w, s]]], a_return)
		end
	end

	#Station to station arcs
	if stationtostation_flag == 1
		for p in pods, w1 in workstations, w2 in [w for w in workstations if w>w1], t in 0:tstep:horizon-arclength[w1,w2]
			a_leave = arcs[nodes[w1, t], nodes[w2, t + arclength[w1, w2]]]
			a_return = arcs[nodes[w2, t], nodes[w1, t + arclength[w2, w1]]]
			
			push!(podarcset[p], a_leave)
			push!(podarcset[p], a_return)
			push!(arcpodset[a_leave], p)
			push!(arcpodset[a_return], p)
			push!(A_plus_p[p, nodes[w1, t]], a_leave)
			push!(A_plus_p[p, nodes[w2, t]], a_return)
			push!(A_minus_p[p, nodes[w2, t + arclength[w1, w2]]], a_leave)
			push!(A_minus_p[p, nodes[w1, t + arclength[w2, w1]]], a_return)
		end
	end

	#Stationary arcs
	if anystoragelocation_flag == 1
		for p in pods, s in storagelocs, t in 0:tstep:horizon-tstep
			a = arcs[nodes[s, t], nodes[s, t + tstep]]
			
			push!(podarcset[p], a)
			push!(arcpodset[a], p)
			push!(A_plus_p[p, nodes[s, t]], a)
			push!(A_minus_p[p,nodes[s, t + tstep]], a)

			for w in workstations
				a = arcs[nodes[w, t], nodes[w, t + tstep]]
				
				push!(podarcset[p], a)
				push!(arcpodset[a], p)
				push!(A_plus_p[p, nodes[w, t]], a)
				push!(A_minus_p[p,nodes[w, t + tstep]], a)
			end
		end
	else
		for p in pods, t in 0:tstep:horizon-tstep
			s = podstorageloc[p]
			a = arcs[nodes[s, t], nodes[s, t + tstep]]
			
			push!(podarcset[p], a)
			push!(arcpodset[a], p)
			push!(A_plus_p[p, nodes[s, t]], a)
			push!(A_minus_p[p,nodes[s, t + tstep]], a)

			for w in workstations
				a = arcs[nodes[w, t], nodes[w, t + tstep]]
				
				push!(podarcset[p], a)
				push!(arcpodset[a], p)
				push!(A_plus_p[p, nodes[w, t]], a)
				push!(A_minus_p[p,nodes[w, t + tstep]], a)
			end
		end
	end

	#Extended arcs
	for p in pods, w in workstations, t in [dummystarttime]
		s = podstorageloc[p]
		podtostationtime = arclength[s, w]
		for t2 in dummystarttime+tstep:tstep:dummystarttime+podtostationtime
			a_leave = extendedarcs[extendednodes[s, t], extendednodes[w, t2]]
			a_return = extendedarcs[extendednodes[w, t], extendednodes[s, t2]]
			
			push!(podarcset[p], a_leave)
			push!(podarcset[p], a_return)
			push!(arcpodset[a_leave], p)
			push!(arcpodset[a_return], p)
			push!(A_plus_p[p, extendednodes[s, t]], a_leave)
			push!(A_plus_p[p, extendednodes[w, t]], a_return)
			push!(A_minus_p[p, extendednodes[w, t2]], a_leave)
			push!(A_minus_p[p, extendednodes[s, t2]], a_return)
		end
	end
	if anystoragelocation_flag == 1
		for p in pods, w in workstations, t in [dummyendtime], s in storagelocs
			podtostationtime = arclength[s, w]
			for t1 in dummyendtime-podtostationtime:tstep:horizon
				a_leave = extendedarcs[extendednodes[s, t1], extendednodes[w, t]]
				a_return = extendedarcs[extendednodes[w, t1], extendednodes[s, t]]
				
				push!(podarcset[p], a_leave)
				push!(podarcset[p], a_return)
				push!(arcpodset[a_leave], p)
				push!(arcpodset[a_return], p)
				push!(A_plus_p[p, extendednodes[s, t1]], a_leave)
				push!(A_plus_p[p, extendednodes[w, t1]], a_return)
				push!(A_minus_p[p, extendednodes[w, t]], a_leave)
				push!(A_minus_p[p, extendednodes[s, t]], a_return)
			end
		end
		for p in pods, t in [dummystarttime, horizon], s in storagelocs
			a = extendedarcs[extendednodes[s, t], extendednodes[s, t + tstep]]
			
			push!(podarcset[p], a)
			push!(arcpodset[a], p)
			push!(A_plus_p[p, extendednodes[s, t]], a)
			push!(A_minus_p[p,extendednodes[s, t + tstep]], a)

			for w in workstations
				a = extendedarcs[extendednodes[w, t], extendednodes[w, t + tstep]]
				
				push!(podarcset[p], a)
				push!(arcpodset[a], p)
				push!(A_plus_p[p, extendednodes[w, t]], a)
				push!(A_minus_p[p,extendednodes[w, t + tstep]], a)
			end
		end
	else
		for p in pods, w in workstations, t in [dummyendtime]
			s = podstorageloc[p]
			podtostationtime = arclength[s, w]
			for t1 in dummyendtime-podtostationtime:tstep:horizon
				a_leave = extendedarcs[extendednodes[s, t1], extendednodes[w, t]]
				a_return = extendedarcs[extendednodes[w, t1], extendednodes[s, t]]
				
				push!(podarcset[p], a_leave)
				push!(podarcset[p], a_return)
				push!(arcpodset[a_leave], p)
				push!(arcpodset[a_return], p)
				push!(A_plus_p[p, extendednodes[s, t1]], a_leave)
				push!(A_plus_p[p, extendednodes[w, t1]], a_return)
				push!(A_minus_p[p, extendednodes[w, t]], a_leave)
				push!(A_minus_p[p, extendednodes[s, t]], a_return)
			end
		end
		for p in pods, t in [dummystarttime, horizon]
			s = podstorageloc[p]
			a = extendedarcs[extendednodes[s, t], extendednodes[s, t + tstep]]
			
			push!(podarcset[p], a)
			push!(arcpodset[a], p)
			push!(A_plus_p[p, extendednodes[s, t]], a)
			push!(A_minus_p[p,extendednodes[s, t + tstep]], a)

			for w in workstations
				a = extendedarcs[extendednodes[w, t], extendednodes[w, t + tstep]]
				
				push!(podarcset[p], a)
				push!(arcpodset[a], p)
				push!(A_plus_p[p, extendednodes[w, t]], a)
				push!(A_minus_p[p,extendednodes[w, t + tstep]], a)
			end
		end
	end
	
	return pods, podarcset, A_minus_p, A_plus_p, podnodeset, arcpodset

end

#-----------------------------------------------------------------------------------#

function getorders(orderheader, orderdetail, allitems, horizonstart, horizonend)

	orders, orders_trueindex, ordersdue = [], [], []
	arrival, deadline = Dict(), Dict()
	itemson = Dict()
	ordermap = Dict()

	orderindex = 1
	for row in 1:size(orderheader)[1]
		#if (orderheader[row,3] <= horizonstart) & (orderheader[row,3] <= horizonend) & (orderheader[row,4] >= horizonstart + 600)
		if (orderheader[row,3] <= horizonstart) & (orderheader[row,3] >= horizonend - horizon*4)
			push!(orders_trueindex, orderheader[row,1])
			ordermap[orderheader[row,1]] = orderindex
			push!(orders, orderindex)

			arrival[orderindex] = tstep * round((orderheader[row,3] - horizonstart) / tstep, digits=0) 
			#deadline[orderheader[row,1]] = tstep * round((orderheader[row,4] - horizonstart) / tstep, digits=0)
			
			#What should this be???????
			deadline[orderindex] = max(tstep * round((orderheader[row,4] - horizonstart) / tstep, digits=0) , rand([240,300,360]))
			
			if tstep * round((orderheader[row,4] - horizonstart) / tstep, digits=0) <= horizon
				push!(ordersdue, orderindex)
			end

			orderindex += 1
		end
	end

	for m in orders
		itemson[m] = []
	end

	for row in 1:size(orderdetail)[1]
		if orderdetail[row,1] in orders_trueindex
			m, i = ordermap[convert(Int64, orderdetail[row,1])], convert(Int64, orderdetail[row,2])
			push!(itemson[m], i)
		end
	end	

	return orders, itemson, deadline, ordersdue		

end

#-----------------------------------------------------------------------------------#

function configureallcapacities(workstationordercapacity, throughputperstation, intersectioncapacity, workstationqueuecapacity)

	C, B, queuecap = Dict(), Dict(), Dict()
	intersectionmaxpods = Dict()
	for w in workstations
		C[w] = workstationordercapacity 
		B[w] = throughputperstation * (tstep/60)
		queuecap[w] = workstationqueuecapacity
	end
	for l in intersections
		intersectionmaxpods[l] = intersectioncapacity
	end

	return C, B, intersectionmaxpods, queuecap

end

#-----------------------------------------------------------------------------------#

function getpodsstoredat()

	podsstoredat = Dict()
	for s in storagelocs
		podsstoredat[s] = []
	end
	for p in allpods
		push!(podsstoredat[podstorageloc[p]], p)
	end

	return podsstoredat

end

#-----------------------------------------------------------------------------------#

function getobjective()

	C=Dict()
	for w in workstations
		C[w] = workstationordercapacity
	end

	return C

end

#-----------------------------------------------------------------------------------#

function calcwarehousedistances()

	#Calculate distance between all locations
	warehousedistance = Dict()
	for l1 in 1:numstoragelocs+numworkstations, l2 in l1:numstoragelocs+numworkstations
		warehousedistance[l1, l2] = abs(loccoords[l1,1] - loccoords[l2,1]) + abs(loccoords[l1,2] - loccoords[l2,2])
		warehousedistance[l2, l1] = warehousedistance[l1, l2]
	end

	return warehousedistance

end

#-----------------------------------------------------------------------------------#

#This is a TEMPORARY FIX - Long term, we need to reconsider how we index variable h_impwt when there are two (or more) of item i on order m

function orderitemtempfix(itemson)
	
	currorderlist = Dict(first(orders) => [])
	for m in orders
		for i in itemson[m]
			try
				if i in currorderlist[m]
					#If there are multiple of one item i on order m, change the second to be another randomly selected item
					newitemnum = convert(Int64,rand(setdiff(items, currorderlist[m])))
					push!(currorderlist[m], newitemnum)
				else
					push!(currorderlist[m], i)
				end
			catch
			    currorderlist[m] = [i]
			end
		end
	end

	for m in orders
		itemson[m] = copy(currorderlist[m])
	end

	orderswith = Dict()
	for i in allitems
		orderswith[i] = []
	end
	for m in orders, i in itemson[m]
		push!(orderswith[i], m)
	end

	items = []
	for i in allitems
		if orderswith[i] != []
			push!(items, i)
		end
	end

	return itemson, items, orderswith

end

#-----------------------------------------------------------------------------------#