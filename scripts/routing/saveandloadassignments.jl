
tstep_cong = 15
tstep_r = 30
maxtraveltime = tstep_r * ceil(((warehouse_x_length_meters + warehouse_y_length_meters) / podspeed) / tstep_r)
horizon_r = horizon*2
routenodelookup, routenodes, routenodes_end, numnodes_r, extendednumnodes_r = createroutingnodes()
routearcs, routearclookup, numarcs_r, RA_plus, RA_minus, newarclength, routearcs_path, routearclookup_path, A_space = createroutingarcs()

function savetaskassignments(assignmentfilename)

	orders_r, pods_r = [], []
	workstationassignment, ordersassignedto, orderassignmentbegan = Dict(), Dict(), Dict()
	podswith_r, podsfor = Dict(), Dict()
	orderopentime, orderclosetime = Dict(), Dict()
	poddelaypenalty = Dict()

	for s in 1:numpartitions
		currpartition, currsol = partitioninfo[s], partitionsolution[s]

		for w in currpartition.workstations
			ordersassignedto[w] = []
		end
		for w in currpartition.workstations, t in times, (m,i,p) in currsol.itempodpicklist[w,t]
			orders_r = union(orders_r, m)
			pods_r = union(pods_r, p)
			workstationassignment[m] = w
			try 
				orderopentime[m] = min(t,orderopentime[m])
			catch
				orderopentime[m] = t
			end
			orderclosetime[m] = horizon_r
			push!(ordersassignedto[w], m)
		end
		for m in orders_r, i in itemson[m]
			podswith_r[i,m] = []
		end
		for m in orders_r
			podsfor[m] = []
		end
		itemsfrom = Dict()
		for w in currpartition.workstations, t in times, (m,i,p) in currsol.itempodpicklist[w,t]
			try
				#println("$w, $t ==> ", itempodpicklist[w,t])
				push!(podswith_r[i,m], p)
				podsfor[m] = union(podsfor[m], p)
				try 
					push!(itemsfrom[p,m], i)
				catch
					itemsfrom[p,m] = [i]
				end
			catch
				println("Removed $m, $i, $p, $w, $t")
			end
		end
		for m in orders_r, t in times
			orderassignmentbegan[m,t] = 0
		end

		for m in orders_r
			for t in times
				if m in currsol.ordersopen[workstationassignment[m], t] 
					for t2 in t:tstep:horizon_r
						orderassignmentbegan[m,t2] = 1
					end
					try 
						orderopentime[m] = min(t,orderopentime[m])
					catch
						orderopentime[m] = t
					end
					for t2 in t:tstep:horizon
						
						if !(m in currsol.ordersopen[workstationassignment[m], t2])
							orderclosetime[m] = t2
							break
						end
					end
					break
				end
			end	
		end

		
		for m in orders_r, p in podsfor[m], t in orderopentime[m]:tstep_r:horizon_r
			n = routenodes[workstationassignment[m], t]
			for a in RA_plus[n]
				poddelaypenalty[m,p,a] = max(0, t - orderclosetime[m]) / tstep
			end
		end
		for a in setdiff(1:numarcs_r, A_space), m in orders_r, p in podsfor[m]
			poddelaypenalty[m,p,a] = 0
		end
	end

	save(assignmentfilename, Dict("orders_r" => orders_r, "pods_r" => pods_r, "podswith_r" => podswith_r, "workstationassignment" => workstationassignment, "ordersassignedto" => ordersassignedto, "orderassignmentbegan" => orderassignmentbegan) )

	return orders_r, pods_r, podswith_r, workstationassignment, ordersassignedto, orderassignmentbegan, orderopentime, orderopentime, podsfor, poddelaypenalty, itemsfrom

end

#-----------------------------------------------------------------------------------#

function loadinstance(instancefilename)
	
	data = load(instancefilename)

	orders_r = data["orders_r"]
	pods_r = data["pods_r"]
	podswith_r = data["podswith_r"]
	workstationassignment = data["workstationassignment"]
	ordersassignedto = data["ordersassignedto"]
	orderassignmentbegan = data["orderassignmentbegan"]

	return orders_r, pods_r, podswith_r, workstationassignment, ordersassignedto, orderassignmentbegan

end

#-----------------------------------------------------------------------------------#

function loadtaskassignments(assignmentfilename)
	
	data = load(assignmentfilename)

	orders_r = data["orders_r"]
	pods_r = data["pods_r"]
	podswith_r = data["podswith_r"]
	workstationassignment = data["workstationassignment"]
	ordersassignedto = data["ordersassignedto"]
	orderassignmentbegan = data["orderassignmentbegan"]

	return orders_r, pods_r, podswith_r, workstationassignment, ordersassignedto, orderassignmentbegan

end
