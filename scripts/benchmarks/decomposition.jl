
function getorderpodlists()
	
	podsfor = Dict()

	for m in orders
		podsfor[m] = []
		for i in itemson[m]
			podsfor[m] = vcat(podsfor[m], podswith[i])
		end
	end

	return podsfor

end

#------------------------------------------------------------------------------------------#

function calculateordercompatibility()
	
	podsfor = getorderpodlists()

	compat = Dict()

	for i1 in 1:length(orders), i2 in i1+1:length(orders)
		m1, m2 = orders[i1], orders[i2]
		compat[m1,m2] = length(intersect(podsfor[m1], podsfor[m2]))
		compat[m2,m1] = compat[m1,m2]
	end

	return compat

end

#------------------------------------------------------------------------------------------#

function calculateorderpodcompatibility()

	compat = Dict()

	for m in orders, p in pods
		compat[m,p] = length(intersect(itemson[m], podstartinventory[p]))
	end

	return compat

end

#------------------------------------------------------------------------------------------#

function getorderworkstationassignment(z)

	ordersassignedto, stationassign = Dict(), Dict()
	for w in workstations
		ordersassignedto[w] = []
		for m in orders
			if getvalue(z[m,w]) > 1e-4
				push!(ordersassignedto[w], m)
				stationassign[m] = w
			end
		end
	end

	return ordersassignedto, stationassign

end

#------------------------------------------------------------------------------------------#

function orderworkstationmodel(ipoutputflag=1, iptimelimit=60) #*10)

	compat = calculateordercompatibility()

	totalitems = sum(length(itemson[m]) for m in orders)
	item_max, item_min = (totalitems/numworkstations)*1.1, (totalitems/numworkstations)*0.9

	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", ipoutputflag)
	set_optimizer_attribute(model, "MIPGap", 0.01)

	#Variables
	@variable(model, z[orders, workstations], Bin)
	@variable(model, ordersync[orders, orders, workstations], Bin)

	#Objective
	@objective(model, Max, sum(sum(sum(compat[m1,m2] * ordersync[m1,m2,w] for m2 in setdiff(orders,m1)) for m1 in orders) for w in workstations))

	#Order-workstation assignment 
	@constraint(model, stationassignment[m = orders], sum(z[m, w] for w in workstations) == 1)
	@constraint(model, ordermaximum[w = workstations], sum(length(itemson[m]) * z[m, w] for m in orders) <= item_max)
	@constraint(model, orderminimum[w = workstations], sum(length(itemson[m]) * z[m, w] for m in orders) >= item_min)
	@constraint(model, ordersynccongif[m1 in orders, m2 in setdiff(orders,m1), w in workstations], ordersync[m1,m2,w] <= 0.5*z[m1,w] + 0.5*z[m2,w])

	#Solve IP
	status = optimize!(model)

	if termination_status(model) == MOI.OPTIMAL
		solvetime = solve_time(model)
		obj = getobjectivevalue(model)
		ordersassignedto, stationassign = getorderworkstationassignment(z)
	elseif termination_status(model) == MOI.TIME_LIMIT
		println("Time out")
		solvetime = solve_time(model)
		obj = getobjectivevalue(model)
		ordersassignedto, stationassign = getorderworkstationassignment(z)
	else
		println("Error in solving!")
		solvetime = solve_time(model)
		obj = 0
	end

	#====================================================#

	return obj, solvetime, z, ordersassignedto, stationassign

end

#------------------------------------------------------------------------------------------#

function getpodassignments(h)

	podassign = Dict()
	for m in orders, i in itemson[m]
		for p in podswith[i]
			if getvalue(h[m,i,p]) > 1e-4
				podassign[i,m] = p
				break
			end
		end
	end

	return podassign

end

#------------------------------------------------------------------------------------------#

#Which pods will I use to fulfill each order? Minimize number of pods used at each workstation 
#Objective here is to maximize the opportunities for pod overlap at the workstations by fulfilling a lot of orders from the same pods. We'll work out the pick/assembly timing details later. 
function orderitempodmodel(ordersassignedto, ipoutputflag=1, iptimelimit=60) #*60*3)

	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", ipoutputflag)
	set_optimizer_attribute(model, "MIPGap", 0.01)

	#Variables
	@variable(model, h[m = orders, i = itemson[m], p = podswith[i]], Bin)
	@variable(model, y[ws = workstations, p = pods], Bin)

	#Objective
	@objective(model, Min, sum(sum(y[ws,p] for p in pods) for ws in workstations) )

	#Item-order-pod assignment
	@constraint(model, orderdelivery[m = orders, i = itemson[m]], sum(h[m, i, p] for p in podswith[i]) == 1)

	#Inventory constraints
	@constraint(model, maxinventory[i = items, p = podswith[i]], sum(h[m,i,p] for m in orderswith[i]) <= inventory[i,p])

	#Workstation-pod assignments
	@constraint(model, podstationassign[ws = workstations, p = pods, m in ordersassignedto[ws], i in intersect(itemson[m], podstartinventory[p])], y[ws,p] >= h[m,i,p] )

	#====================================================#

	#Solve IP
	status = optimize!(model)

	if termination_status(model) == MOI.OPTIMAL
		solvetime = solve_time(model)
		obj = getobjectivevalue(model)
		podassign = getpodassignments(h)
	elseif termination_status(model) == MOI.TIME_LIMIT
		println("Time out in the subproblem")
		solvetime = solve_time(model)
		obj = getobjectivevalue(model)
		podassign = getpodassignments(h)
	else
		println("Error in solving!")
		solvetime = solve_time(model)
		obj = 0
	end

	#====================================================#

	return obj, solvetime, h, podassign

end

#------------------------------------------------------------------------------------------#

function podflowmodel(stationassign, podassign, ipoutputflag=1, iptimelimit=60) #*60*24)

	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", ipoutputflag)
	set_optimizer_attribute(model, "MIPGap", 0.01)

	#Variables
	@variable(model, h[m = orders, i = itemson[m], times], Bin)
	@variable(model, y[p = pods, a = podarcset[p]], Bin)
	@variable(model, v[orders, times], Bin)
	@variable(model, f[orders,  times], Bin)
	@variable(model, g[orders, times], Bin)

	#Objective
	@objective(model, Max, sum(sum(h[m, i, horizon] for i in itemson[m]) for m in orders) )

	#Linking constraints																														
	@constraint(model, podlink[m = orders, i = itemson[m], t = setdiff(times,0)], h[m,i,t] - h[m,i,t-tstep] <= sum(y[podassign[i,m],a] for a in setdiff(A_plus_p[podassign[i,m],nodes[stationassign[m],t]], A_queues))) 
	@constraint(model, podlinkzero[m = orders, i = itemson[m]], h[m,i,0] <= sum(y[podassign[i,m],a] for a in setdiff(A_plus_p[podassign[i,m],nodes[stationassign[m],0]], A_queues))) 
	@constraint(model, nondecreasing[m = orders, i = itemson[m], t = setdiff(times,0)], h[m,i,t] >= h[m,i,t-tstep] )

	#Pod movement constraints
	@constraint(model, podstartlocations[p = pods], sum(y[p,a] for a in A_plus_p[p, podstartnode[p]]) == 1)
	@constraint(model, podflowbalance[p = pods, n in setdiff(1:numnodes, union(podstartnode[p], N_end))], sum(y[p,a] for a in A_minus_p[p,n]) - sum(y[p,a] for a in A_plus_p[p,n]) == 0)

	#Workstation capacity constraints
	@constraint(model, fconfig[m = orders, i = itemson[m], t = times], f[m,t] >= h[m,i,t])
	@constraint(model, gconfig[m = orders, i = itemson[m], t = times], g[m,t] <= h[m,i,t])
	@constraint(model, vconfig[m = orders, t = times], v[m,t] == f[m,t] - g[m,t])
	@constraint(model, stationcapacity[w = workstations, t = times], sum(v[m,t] for m in ordersassignedto[w]) <= C[w] )

	#Workstation throughput constraints
	@constraint(model, stationthroughput[w = workstations, t = setdiff(times,0)], itemprocesstime * sum(sum(h[m,i,t] - h[m,i,t-tstep] for i in itemson[m]) for m in ordersassignedto[w]) + podprocesstime * sum(sum(y[p,a] for a in intersect(A_plus_p[p,nodes[w,t]], A_space)) for p in pods) <= tstep)
	@constraint(model, stationthroughput_first[w = workstations], itemprocesstime * sum(sum(h[m,i,0] for i in itemson[m]) for m in ordersassignedto[w]) + podprocesstime * sum(sum(y[p,a] for a in intersect(A_plus_p[p,nodes[w,0]], A_space)) for p in pods) <= tstep)

	#Congestion
	constraint(model, maxcongestion[l in intersections, t in 0:congestiontstep:horizon], sum(sum(congestioncontribution[a,l,t] * y[p,a] for a in intersect(congestionarcs[l,t], podarcset[p])) for p in pods) <= intersectionmaxpods[l] )

	#====================================================#

	#Solve IP
	status = optimize!(model)

	if termination_status(model) == MOI.OPTIMAL
		solvetime = solve_time(model)
		obj = getobjectivevalue(model)
		println("Total items picked = ", sum(sum(getvalue(h[m, i, horizon]) for i in itemson[m]) for m in orders) )
	elseif termination_status(model) == MOI.TIME_LIMIT
		println("Time out in the subproblem")
		solvetime = solve_time(model)
		obj = getobjectivevalue(model)
		println("Total items picked = ", sum(sum(getvalue(h[m, i, horizon]) for i in itemson[m]) for m in orders) )
	else
		println("Error in solving!")
		solvetime = solve_time(model)
		obj = 10000000
		println("Total items picked = ", sum(sum(getvalue(h[m, i, horizon]) for i in itemson[m]) for m in orders) )
	end

	#====================================================#

	#=
	podsegments = []
	for p in pods, a in intersect(union(A_queues, A_space), podarcset[p])
		if getvalue(y[p,a]) > 0.001
			loc1, t1 = nodelookup[arclookup[a][1]]
			loc2, t2 = nodelookup[arclookup[a][2]]
			push!(podsegments, (p, t1, t2, loc1, loc2))
		end
	end

	openorders = Dict()
	for m in orders, w in workstations, t in times
		if (getvalue(v[m,w,t]) > 0.01) && (sum(sum(getvalue(h[m,i,p,w,t]) for p in podswith[i]) for i in itemson[m]) >= 1) && ( sum(sum(getvalue(h[m,i,p,w,t]) for p in podswith[i]) for i in itemson[m]) < length(itemson[m]))
			try
				push!(openorders[nodes[w,t]], m)
			catch
				openorders[nodes[w,t]] = [m]
			end
		end
	end
	=#

	#====================================================#

	return obj, solvetime, h, y, v

end

#------------------------------------------------------------------------------------------#

#Which pods will I use to fulfill each order? Minimize number of pods used at each workstation 
#Objective here is to maximize the opportunities for pod overlap at the workstations by fulfilling a lot of orders from the same pods. We'll work out the pick/assembly timing details later. 
function orderworkstationitempodmodel(ipoutputflag=1, iptimelimit=60) #*60*6)

	compat = calculateordercompatibility()

	totalitems = sum(length(itemson[m]) for m in orders)
	item_max, item_min = (totalitems/numworkstations)*1.1, (totalitems/numworkstations)*0.9

	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", ipoutputflag)
	set_optimizer_attribute(model, "MIPGap", 0.01)

	#Variables
	@variable(model, z[orders, workstations], Bin)
	@variable(model, h[m = orders, i = itemson[m], p = podswith[i]], Bin)
	@variable(model, y[ws = workstations, p = pods], Bin)

	#Objective
	@objective(model, Min, sum(sum(y[ws,p] for p in pods) for ws in workstations) )

	#Order-workstation assignment 
	@constraint(model, stationassignment[m = orders], sum(z[m, w] for w in workstations) == 1)
	@constraint(model, ordermaximum[w = workstations], sum(length(itemson[m]) * z[m, w] for m in orders) <= item_max)
	@constraint(model, orderminimum[w = workstations], sum(length(itemson[m]) * z[m, w] for m in orders) >= item_min)

	#Item-order-pod assignment
	@constraint(model, orderdelivery[m = orders, i = itemson[m]], sum(h[m, i, p] for p in podswith[i]) == 1)

	#Inventory constraints
	@constraint(model, maxinventory[i = items, p = podswith[i]], sum(h[m,i,p] for m in orderswith[i]) <= inventory[i,p])

	#Linking
	@constraint(model, podstationassign[ws in workstations, p in pods, m in orders, i in intersect(itemson[m], podstartinventory[p])], y[ws,p] >= h[m,i,p] + z[m,ws] - 1 )

	#====================================================#

	#Solve IP
	status = optimize!(model)

	if termination_status(model) == MOI.OPTIMAL
		solvetime = solve_time(model)
		obj = getobjectivevalue(model)
		ordersassignedto, stationassign = getorderworkstationassignment(z)
		podassign = getpodassignments(h, z)
	elseif termination_status(model) == MOI.TIME_LIMIT
		println("Time out in the subproblem")
		solvetime = solve_time(model)
		obj = getobjectivevalue(model)
		ordersassignedto, stationassign = getorderworkstationassignment(z)
		podassign = getpodassignments(h)
	else
		println("Error in solving!")
		solvetime = solve_time(model)
		obj = 0
	end

	#====================================================#

	return obj, solvetime, h, z, y, ordersassignedto, stationassign, podassign

end

#------------------------------------------------------------------------------------------#

function getpickschedule(h)

	pickschedule = Dict()

	for m in orders, i in itemson[m]
		pickschedule[m,i,0] = getvalue(h[m,i,0])
	end
	for m in orders, i in itemson[m], t in setdiff(times,0)
		pickschedule[m,i,t] = getvalue(h[m,i,t]) - getvalue(h[m,i,t-tstep])
	end

	return pickschedule

end

#------------------------------------------------------------------------------------------#

function workstationschedulingmodel(stationassign, podassign, ipoutputflag=1, iptimelimit=60) #*60*24)

	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", ipoutputflag)
	set_optimizer_attribute(model, "MIPGap", 0.01)

	#Variables
	@variable(model, h[m = orders, i = itemson[m], times], Bin)
	@variable(model, y[p = pods, a = podarcset[p]], Bin)
	@variable(model, v[orders, times], Bin)
	@variable(model, f[orders,  times], Bin)
	@variable(model, g[orders, times], Bin)

	#Objective
	@objective(model, Max, sum(sum(h[m, i, horizon] for i in itemson[m]) for m in orders) )

	#Linking constraints																														
	@constraint(model, podlink[m = orders, i = itemson[m], t = setdiff(times,0)], h[m,i,t] - h[m,i,t-tstep] <= sum(y[podassign[i,m],a] for a in setdiff(A_plus_p[podassign[i,m],nodes[stationassign[m],t]], A_queues))) 
	@constraint(model, podlinkzero[m = orders, i = itemson[m]], h[m,i,0] <= sum(y[podassign[i,m],a] for a in setdiff(A_plus_p[podassign[i,m],nodes[stationassign[m],0]], A_queues))) 
	@constraint(model, nondecreasing[m = orders, i = itemson[m], t = setdiff(times,0)], h[m,i,t] >= h[m,i,t-tstep] )

	#Pod movement constraints
	@constraint(model, podstartlocations[p = pods], sum(y[p,a] for a in A_plus_p[p, podstartnode[p]]) == 1)
	@constraint(model, podflowbalance[p = pods, n in setdiff(1:numnodes, union(podstartnode[p], N_end))], sum(y[p,a] for a in A_minus_p[p,n]) - sum(y[p,a] for a in A_plus_p[p,n]) == 0)

	#Workstation capacity constraints
	@constraint(model, fconfig[m = orders, i = itemson[m], t = times], f[m,t] >= h[m,i,t])
	@constraint(model, gconfig[m = orders, i = itemson[m], t = times], g[m,t] <= h[m,i,t])
	@constraint(model, vconfig[m = orders, t = times], v[m,t] == f[m,t] - g[m,t])
	@constraint(model, stationcapacity[w = workstations, t = times], sum(v[m,t] for m in ordersassignedto[w]) <= C[w] )

	#Workstation throughput constraints
	@constraint(model, stationthroughput[w = workstations, t = setdiff(times,0)], itemprocesstime * sum(sum(h[m,i,t] - h[m,i,t-tstep] for i in itemson[m]) for m in ordersassignedto[w]) + podprocesstime * sum(sum(y[p,a] for a in intersect(A_plus_p[p,nodes[w,t]], A_space)) for p in pods) <= tstep)
	@constraint(model, stationthroughput_first[w = workstations], itemprocesstime * sum(sum(h[m,i,0] for i in itemson[m]) for m in ordersassignedto[w]) + podprocesstime * sum(sum(y[p,a] for a in intersect(A_plus_p[p,nodes[w,0]], A_space)) for p in pods) <= tstep)

	#====================================================#

	#Solve IP
	status = optimize!(model)

	if termination_status(model) == MOI.OPTIMAL
		solvetime = solve_time(model)
		obj = getobjectivevalue(model)
		pickschedule = getpickschedule(h)
		println("Total items picked = ", sum(sum(getvalue(h[m, i, horizon]) for i in itemson[m]) for m in orders) )
	elseif termination_status(model) == MOI.TIME_LIMIT
		println("Time out in the subproblem")
		solvetime = solve_time(model)
		obj = getobjectivevalue(model)
		pickschedule = getpickschedule(h)
		println("Total items picked = ", sum(sum(getvalue(h[m, i, horizon]) for i in itemson[m]) for m in orders) )
	else
		println("Error in solving!")
		solvetime = solve_time(model)
		obj = 10000000
		println("Total items picked = ", sum(sum(getvalue(h[m, i, horizon]) for i in itemson[m]) for m in orders) )
	end

	#====================================================#

	return obj, solvetime, h, v, pickschedule

end

#------------------------------------------------------------------------------------------#

function congestioncheckmodel(stationassign, podassign, pickschedule, ipoutputflag=1, iptimelimit=60) #*60*24)

	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", ipoutputflag)
	set_optimizer_attribute(model, "MIPGap", 0.01)

	#Variables
	@variable(model, z[m = orders], Bin)
	@variable(model, h[m = orders, i = itemson[m], times], Bin)
	@variable(model, y[p = pods, a = podarcset[p]], Bin)

	#Objective
	@objective(model, Max, sum(sum(h[m, i, horizon] for i in itemson[m]) for m in orders) )

	#Pick assignments
	@constraint(model, picktime[m = orders, i = itemson[m], t = times], h[m,i,t] <= pickschedule[m,i,t])
	@constraint(model, pickwholeorder[m = orders, i = itemson[m]], sum(h[m,i,t] for t in times) == z[m] )

	#Linking constraints																														
	@constraint(model, podlink[m = orders, i = itemson[m], t = times], h[m,i,t] <= sum(y[podassign[i,m],a] for a in setdiff(A_plus_p[podassign[i,m],nodes[stationassign[m],t]], A_queues))) 

	#Pod movement constraints
	@constraint(model, podstartlocations[p = pods], sum(y[p,a] for a in A_plus_p[p, podstartnode[p]]) == 1)
	@constraint(model, podflowbalance[p = pods, n in setdiff(1:numnodes, union(podstartnode[p], N_end))], sum(y[p,a] for a in A_minus_p[p,n]) - sum(y[p,a] for a in A_plus_p[p,n]) == 0)

	#Congestion
	#@constraint(model, maxcongestion[l in intersections, t in 0:congestiontstep:horizon], sum(sum(congestioncontribution[a,l,t] * y[p,a] for a in intersect(congestionarcs[l,t], podarcset[p])) for p in pods) <= intersectionmaxpods[l] )

	#====================================================#

	#Solve IP
	status = optimize!(model)

	if termination_status(model) == MOI.OPTIMAL
		solvetime = solve_time(model)
		obj = getobjectivevalue(model)
		println("Total items picked = ", sum(sum(getvalue(h[m, i, horizon]) for i in itemson[m]) for m in orders) )
	elseif termination_status(model) == MOI.TIME_LIMIT
		println("Time out in the subproblem")
		solvetime = solve_time(model)
		obj = getobjectivevalue(model)
		println("Total items picked = ", sum(sum(getvalue(h[m, i, horizon]) for i in itemson[m]) for m in orders) )
	else
		println("Error in solving!")
		solvetime = solve_time(model)
		obj = 10000000
		println("Total items picked = ", sum(sum(getvalue(h[m, i, horizon]) for i in itemson[m]) for m in orders) )
	end

	#====================================================#

	return obj, solvetime, h, y, z

end

#------------------------------------------------------------------------------------------#

function writedecompositionoutputs(outputfilename, obj_full, solvetime_ow, solvetime_oip, solvetime_t, solvetime_full, h_full, y_full, v_full)

	itemspicked = 0
	itemspickedfrom = Dict()
	podpicklist = Dict()
	for m in orders
		itemspickedfrom[m] = []
	end
	for w in workstations, t in times
		podpicklist[w,t] = []
	end
	for m in orders, i in itemson[m]
		for t in times
			if getvalue(h_full[m,i,t]) > 1e-4
				itemspicked += 1
				push!(itemspickedfrom[m], i)
				push!(podpicklist[stationassign[m],t], podassign[i,m])
			end
		end
	end
	podsused = sum(sum(length(podpicklist[w,t]) for w in workstations) for t in times)

	ordersworked, orderscompleted = [], []
	for m in orders
		if length(itemspickedfrom[m]) == length(itemson[m])
			push!(orderscompleted, m)
			push!(ordersworked, m)
		elseif length(itemspickedfrom[m]) >= 1
			push!(ordersworked, m)
		end
	end

	orderfreq = zeros(11)
	for m in orderscompleted
		l = convert(Int64, length(itemson[m]))
		orderfreq[l] += 1
	end

	order_open_time_per_item = sum(sum(getvalue(v_full[m,t]) for t in times) for m in orderscompleted) / sum(length(itemson[m]) for m in orderscompleted)

	congestionat = Dict()
	for l in intersections, t in 0:congestiontstep:horizon
		congestionat[l,t] = 0
		for p in pods, a in intersect(congestionarcs[l,t], podarcset[p])
			congestionat[l,t] += congestioncontribution[a,l,t] * getvalue(y_full[p,a])
		end	
		congestionat[l,t] = congestionat[l,t] / intersectionmaxpods[l]
	end
	listofcongestions = values(congestionat)

	totaldist, podtrips = 0, 0
	for w in workstations, t in times, p in podpicklist[w,t]
		totaldist += warehousedistance[podstorageloc[p],w]
		podtrips += 1
	end

	df = DataFrame(run_id = [run_id],
			test_instance_id = [instance_id],
			warehouse_id = [warehouse_id], 
			random_seed = [random_seed], 
			decomposition_type = [decompositiontype], 
			objective = [obj_full],
			solve_time = [solvetime_ow + solvetime_oip + solvetime_t + solvetime_full],
			decomp1_time = [solvetime_ow],
			decomp2_time = [solvetime_oip],
			decomp3_time = [solvetime_t],
			decomp4_time = [solvetime_full],
			time_utilization = [(itemspicked*itemprocesstime + podsused*podprocesstime) / (tstep * numworkstations * length(times))],
			throughput_utilization = [(itemspicked*(itemprocesstime+podprocesstime)) / (tstep * numworkstations * length(times))],
			bestthroughput_utilization = [(itemspicked*itemprocesstime) / ((tstep - podprocesstime) * numworkstations * length(times))],
			total_orders_worked = [length(ordersworked)],
			total_orders_completed = [length(orderscompleted)],
			congestion_utilization = [mean(listofcongestions)],
			max_congestion = [maximum(listofcongestions)],
			pods_used = [podsused],
			items_picked_per_pod = [itemspicked / podsused],
			pod_distance_travelled = [totaldist / podtrips],
			order_open_time_per_item = [order_open_time_per_item],
			orders_size_1 = [orderfreq[1]],
			orders_size_2 = [orderfreq[2]],
			orders_size_3 = [orderfreq[3]],
			orders_size_4 = [orderfreq[4]],
			orders_size_5 = [orderfreq[5]],
			orders_size_6 = [orderfreq[6]],
			orders_size_7 = [orderfreq[7]],
			orders_size_8 = [orderfreq[8]],
			orders_size_9 = [orderfreq[9]],
			orders_size_10 = [orderfreq[10]],
			orders_size_large = [orderfreq[11]]
           )

	CSV.write(outputfilename, df)

end

#------------------------------------------------------------------------------------------#
