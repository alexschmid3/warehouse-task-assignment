
function createemptysolution(orders, pods, podswith, workstations)

	y_currsol = Dict()
	h_currsol = Dict()
	y_currpath = Dict()
	v_currsol = Dict()

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

	congestionat = Dict()

	solution = (y=y_currsol, h=h_currsol, v=v_currsol, ypath=y_currpath, podsworkedat=podsworkedat, ordersopen=ordersopen, unassignedorders=unassignedorders, itempodpicklist=itempodpicklist, congestionat=congestionat)

	return solution

end