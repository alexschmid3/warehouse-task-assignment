
function updatecongestionat(y_currsol)

	for l in reducedintersections, t in 0:congestiontstep:horizon
		if congestionarcs[l,t] != []
			congestionat[l,t] = sum(sum(congestioncontribution[a,l,t] * y_currsol[p,a] for p in arcpodset[a]) for a in congestionarcs[l,t])
		else 
			congestionat[l,t] = 0
		end
	end

end

#-----------------------------------------------------------------------------#

function update_h(sp, currsol, h)

	#Update order-item-pod-workstation-time assignment
	for m in sp.orders, i in sp.itemson[m], p in sp.podswith[i], w in sp.workstations
		for t in sp.times
			if value(h[m, i, p, w, t]) > 1e-2
				currsol.h[m, i, p, w, t] = value(h[m, i, p, w, t])
				for t2 in t:tstep:horizon
					currsol.h[m, i, p, w, t2] = value(h[m, i, p, w, t])
				end
				currsol.stationassign[m] = w
				df = DataFrame(iteration=[counter], order=[m], station=[w])
				break
			elseif (value(h[m, i, p, w, t]) <= 1e-2) & (currsol.h[m, i, p, w, t] >= 1e-2)
				currsol.h[m, i, p, w, t] = value(h[m, i, p, w, t])
				for t2 in t:tstep:horizon
					currsol.h[m, i, p, w, t2] = 0
				end
				#currsol.stationassign[m] = 0
				#break
			else
				currsol.h[m, i, p, w, t] = value(h[m, i, p, w, t])
			end
		end
	end

	return currsol

end

#-----------------------------------------------------------------------------#

function update_y(sp, currsol, y)

	podpath = Dict()
	for p in sp.pods
		podpath[p] = []
	end

	for p in sp.pods, a in sp.podarcset[p]
		currsol.y[p, a] = value(y[p, a])
		if value(y[p, a]) > 0.01
			push!(currsol.ypath[p], a)
			if a <= numarcs
				push!(podpath[p], a)
			end
		else
			remove!(currsol.ypath[p], a)
		end
	end

	podbusy = Dict()
	for p in sp.pods
		podbusy[p] = []
		for a in currsol.ypath[p]
			startloc, endloc = nodelookup[arclookup[a][1]][1], nodelookup[arclookup[a][2]][1]
			starttime, endtime = nodelookup[arclookup[a][1]][2], nodelookup[arclookup[a][2]][2]
			#for w in unique!(setdiff([startloc, endloc], podstorageloc[p])), t in starttime:tstep:endtime
			for w in unique!(intersect([startloc, endloc], workstations)), t in starttime:tstep:endtime
				push!(podbusy[p], (w,t))
			end
		end
		currsol.podbusy[p] = unique!(podbusy[p])
	end

	return currsol, podpath

end

#-----------------------------------------------------------------------------#

function update_v(sp, partition, currsol, v)

	#Update open orders
	for w in sp.workstations, t in sp.times 
		currsol.ordersopen[w,t] = setdiff(currsol.ordersopen[w,t], sp.orders)
	end

	#Need to change to detailed accounting for v, in case an order is marked open when it is complete because capacity constraints aren't binding
	#for m in sp.orders, w in sp.workstations, t in sp.times 
	#	currsol.v[m, w, t] = value(v[m, w, t])
	#	if value(v[m, w, t]) > 1e-4			
	#		push!(currsol.ordersopen[w,t], m)
	#	end
	#end

	#Need to change to detailed accounting for v, in case an order is marked open when it is complete because capacity constraints aren't binding
	for m in sp.orders, w in sp.workstations
		#Use this version for julia 1.6 and higher:
		#if sum(sum(sum(currsol.h[m,i,p,w,t] for t in times) for p in partition.podswith[i]; init=0) for i in itemson[m]; init=0) > 1e-4
		#Use this version for julia 1.5 and lower:
		mysum = 0
		for i in itemson[m], p in partition.podswith[i], t in times
			mysum += currsol.h[m,i,p,w,t]
		end
		if mysum > 1e-4
			mint, maxt = horizon*2, -1
			for i in itemson[m]
				if sum(sum(currsol.h[m,i,p,w,t] for p in partition.podswith[i]) for t in times) < 1e-4
					#Item never delivered, so the order is still open at the end of the time horizon
					maxt = max(maxt, horizon + tstep)
				else
					#Item delivered, so we need to find the delivery time and adjust the min and max delivery times accordingly
					for t in times
						if sum(currsol.h[m,i,p,w,t] for p in partition.podswith[i]) > 1e-4
							mint = min(mint, t)
							maxt = max(maxt, t)
							break
						end
					end
				end
			end
			for t in mint:tstep:maxt - tstep
				currsol.v[m,w,t] = 1
				push!(currsol.ordersopen[w,t], m)
			end
		else
			for t in sp.times
				currsol.v[m,w,t] = 0
			end
		end
	end

	return currsol

end

#-----------------------------------------------------------------------------#

function update_picklists(sp, currsol, h)

	#Update item-pod pick list (used for inventory calculation)
	truepicks = value.(h)
	for w in sp.workstations, t in sp.times
		newlist = []
		for (m,i,p) in currsol.itempodpicklist[w,t]
			if !(m in sp.orders) || !(p in sp.pods) 
				push!(newlist, (m,i,p))
			end
		end	
		currsol.itempodpicklist[w,t] = newlist
		currsol.podsworkedat[w,t] = setdiff(currsol.podsworkedat[w,t], sp.pods)
	end

	for w in sp.workstations, m in sp.orders, i in sp.itemson[m], p in sp.podswith[i]
		if truepicks[m,i,p,w,sp.times[1]] == 1
			push!(currsol.itempodpicklist[w,sp.times[1]], (m, i, p))
			currsol.podsworkedat[w,sp.times[1]] = union(currsol.podsworkedat[w,sp.times[1]], p)
		end
	end

	for w in sp.workstations, m in sp.orders, i in sp.itemson[m], p in sp.podswith[i], t in setdiff(sp.times,sp.times[1])
		if truepicks[m,i,p,w,t] - truepicks[m,i,p,w,t-tstep] == 1
			push!(currsol.itempodpicklist[w,t], (m, i, p))
			currsol.podsworkedat[w,t] = union(currsol.podsworkedat[w,t], p)
		end
	end

	return currsol

end

#-----------------------------------------------------------------------------#

function updateorders(currpartition, currsol)

	#Update unassigned orders
	for m in currpartition.orders
		remove!(currsol.unassignedorders, m) 			#Clear all orders
		push!(currsol.unassignedorders, m)   			#Reset with all orders
	end

	for w in currpartition.workstations, t in times, m in currsol.ordersopen[w,t]
		remove!(currsol.unassignedorders, m)			#Remove the orders being worked
	end

	return currsol

end

#-----------------------------------------------------------------------------#

function updatecongestion(sp, currpartition, podpath)

	#Identify intersections that could be impacted by the subproblem re-optimization
	cong_intersections = [maps.mapintersectiontorow[i] for i in currpartition.intersections]

	for p in sp.pods
		
		#Find the times the current pod was impacted by our re-optimization
		podtimes = []
		for a in podpath[p]
			t1, t2 = max(0, nodelookup[arclookup[a][1]][2]), min(horizon, nodelookup[arclookup[a][2]][2])
			for t in t1:congestiontstep:t2
				push!(podtimes, t)
			end
		end
		podtimes = unique(podtimes)
		cong_times = [maps.maptimetocolumn[t] for t in podtimes]

		#Reset the congestion of that pod for the impacted times - unnecessary
		#currcong[p][cong_times, cong_intersections] = zeros(length(cong_times), length(cong_intersections))

		#Add the congestion generated by the new pod path
		if length(podpath[p]) > 0
			currcong[p][cong_intersections, cong_times] = sum(congestionsignature[a][cong_intersections, cong_times] for a in podpath[p])
		else
			currcong[p][cong_intersections, cong_times] = zeros(length(cong_intersections), length(cong_times))
		end

	end

end

#-----------------------------------------------------------------------------#

function updatesolution(sp, currsol, currpartition, new_obj, h, y, v)

	#Calculate old objective
	old_obj = 0
	for w in sp.workstations, m in sp.orders, i in sp.itemson[m], p in sp.podswith[i] 
		old_obj += currsol.h[m, i, p, w, last(sp.times)] 
	end

	#Update 
	if new_obj - old_obj > 1e-4
		currsol = update_h(sp, currsol, h)
		currsol, podpath = update_y(sp, currsol, y)
		currsol = update_v(sp, currpartition, currsol, v)
		currsol = update_picklists(sp, currsol, h)
		currsol = updateorders(currpartition, currsol)
		updatecongestion(sp, currpartition, podpath)
		println("Solution improved from $old_obj to $new_obj")
	else
		println("Solution remains at $old_obj")
	end

	return currsol

end

# sp, currsol, currpartition, new_obj, h, y, v = sp, currsol, currpartition, sp_obj, h_sp, y_sp, v_sp




