
function createemptysolution(orders, pods, podswith, workstations)

	y_currsol = Dict()
	h_currsol = Dict()
	y_currpath = Dict()
	v_currsol = Dict()
	stationassign = Dict()

	for p in pods, a in podarcset[p]
		y_currsol[p,a] = 0
	end

	for m in orders, i in itemson[m], p in podswith[i], w in workstations, t in extendedtimes
		h_currsol[m, i, p, w, t] = 0
	end
	
	for m in orders, w in workstations, t in extendedtimes
		v_currsol[m,w,t] = 0
	end

	for p in pods
		y_currpath[p] = []
		for t in dummystarttime:tstep:horizon
			a = extendedarcs[extendednodes[podstorageloc[p], t], extendednodes[podstorageloc[p], t+tstep]]
			push!(y_currpath[p], a)
			y_currsol[p,a] += 1
		end
	end

	podsworkedat = Dict()
	itempodpicklist = Dict()
	ordersopen = Dict()
	unassignedorders = copy(orders)

	for w in workstations, t in times
		podsworkedat[w,t] = []
	end

	for w in workstations, t in times
		ordersopen[w,t] = []
	end

	for w in workstations, t in times
		itempodpicklist[w,t] = []
	end

	for m in orders
		stationassign[m] = 0
	end

	podbusy = Dict()
	for p in pods
		podbusy[p] = []
		for a in y_currpath[p]
			startloc, endloc = nodelookup[arclookup[a][1]][1], nodelookup[arclookup[a][2]][1]
			starttime, endtime = nodelookup[arclookup[a][1]][2], nodelookup[arclookup[a][2]][2]
			for w in unique!(setdiff([startloc, endloc], podstorageloc[p])), t in starttime:tstep:endtime
				push!(podbusy[p], (w,t))
			end
		end
		podbusy[p] = unique!(podbusy[p])
	end

	solution = (y=y_currsol, h=h_currsol, v=v_currsol, ypath=y_currpath, podsworkedat=podsworkedat, ordersopen=ordersopen, unassignedorders=unassignedorders, itempodpicklist=itempodpicklist, stationassign=stationassign, podbusy=podbusy)

	return solution

end