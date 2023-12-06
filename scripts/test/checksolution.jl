
function checksolution(currpartition, currsol, printstatements)

	if printstatements == 1
		for w in currpartition.workstations, t in times
			println("---- $w, $t ----")
			for (m,i,p) in currsol.itempodpicklist[w,t]
				println("($m, $i, $p)")
			end
			println("Pods = ", currsol.podsworkedat[w,t])
		end
	end

	#Check whether the throughput constraint is violated and the pods are accounted for
	for w in currpartition.workstations, t in times
		podlist = []
		for (m,i,p) in currsol.itempodpicklist[w,t]
			podlist = union(podlist, p)
		end
		usedtime = length(currsol.itempodpicklist[w,t]) * itemprocesstime + length(podlist) * podprocesstime

		if printstatements == 1
			println("Check 0 - $w, $t --> ", usedtime, " from ", length(currsol.itempodpicklist[w,t]), " items and ", length(podlist), " pods")
		end
		@assert usedtime <= tstep
		if printstatements == 1
			println("Check 0.1 - $w, $t --> ", podlist, " = ", currsol.podsworkedat[w,t])
		end
		@assert podlist == currsol.podsworkedat[w,t]
	end
	
	#Cross-reference h with item pod pick list
	timesassigned = Dict()
	for m in currpartition.orders, i in itemson[m]
		timesassigned[m,i] = []
	end
	for w in currpartition.workstations, t in times
		for (m,i,p) in currsol.itempodpicklist[w,t]

			push!(timesassigned[m,i], (p,w,t))
			
			if printstatements == 1
				println("Check 1 - $m, $i, $p, $w, $t")
			end
			@assert currsol.h[m,i,p,w,t] >= 1e-4
			
			if t >= 0
				if printstatements == 1
					println("Check 1.1 - $m, $i, $p, $w, ", t-tstep)
				end
				@assert currsol.h[m,i,p,w,t-tstep] == 0
			end
		end
	end

	#Ensure no item assigned twice
	for m in currpartition.orders, i in itemson[m]
		if printstatements == 1
			println("Check 1.2 - $m, $i --> ", timesassigned[m,i])
		end
		@assert length(timesassigned[m,i]) <= 1
	end

	#Reverse cross check of picklist with h
	for m in currpartition.orders, i in itemson[m], p in currpartition.podswith[i], w in currpartition.workstations
		for t in times
			if currsol.h[m,i,p,w,t] >= 1e-4
				if printstatements == 1
					println("Check 1.3 - $m, $i, $p, $w, $t")
				end
				@assert (m,i,p) in currsol.itempodpicklist[w,t]
				break
			end
		end
	end

	#Cross-reference v with ordersopen
	for w in currpartition.workstations, t in times
		for m in currsol.ordersopen[w,t]
			if printstatements == 1
				println("Check 2 - $m, $w, $t")
			end
			@assert currsol.v[m,w,t] == 1
		end
	end

	#Cross-reference h with v and station assignment
	for m in currpartition.orders
		numassignedstations = 0
		for w in currpartition.workstations
			if sum(sum(sum(currsol.h[m,i,p,w,t] for t in times) for p in currpartition.podswith[i]) for i in itemson[m]) > 1e-4
				
				#Check stored station assignment
				if printstatements == 1
					println("Check 3 - currsol.stationassign[$m] = ", currsol.stationassign[m], " == $w")
				end
				@assert currsol.stationassign[m] == w
				numassignedstations += 1

				#Check order open assignment
				mint, maxt = horizon*2, -1
				for i in itemson[m]
					for t in times
						if sum(currsol.h[m,i,p,w,t] for p in currpartition.podswith[i]) > 1e-4
							mint = min(mint, t)
							maxt = max(maxt, t)
							break
						end
					end
				end
				for t in mint:tstep:maxt - tstep
					if printstatements == 1
						println("Check 4 - $m, $w, $t")
					end
					@assert currsol.v[m,w,t] == 1
				end

			end
		end
		if printstatements == 1
			println("Check 5 - $m")
		end
		@assert numassignedstations <= 1
	end

	#Check that all items from every order are delivered or the order is left open at the end of the window
	for m in currpartition.orders
		deliveredlist = []
		for i in itemson[m]
			push!(deliveredlist, sum(sum(currsol.h[m,i,p,w,horizon] for p in currpartition.podswith[i]) for w in currpartition.workstations))
		end
		stillopenflag = sum(currsol.v[m,w,horizon] for w in currpartition.workstations)
		fulldeliveredflag = minimum(deliveredlist)
		partiallydeliveredflag = maximum(deliveredlist)

		deliveryokay = 0
		if fulldeliveredflag == 1
			deliveryokay += 1
		elseif (partiallydeliveredflag == 1) & (stillopenflag == 1)
			deliveryokay += 1
		elseif (partiallydeliveredflag == 0) & (stillopenflag == 0)
			deliveryokay += 1
		end

		if printstatements == 1
			println("Check 6 - $m ($fulldeliveredflag, $partiallydeliveredflag, $stillopenflag)")
		end
		@assert deliveryokay == 1
	end
	
	remainingcongestionspace = intersectiontimemaxpods - sum(currcong[p] for p in pods)
	for l in 1:length(intersections), t in 1:length(0:congestiontstep:horizon)
		println("Check 7 - $l, $t - ", remainingcongestionspace[l,t])
		@assert remainingcongestionspace[l,t] >= -1e-4
	end

	println("Solution check passed")

end