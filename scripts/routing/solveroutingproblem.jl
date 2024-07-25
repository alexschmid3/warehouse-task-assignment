
function solveroutingproblem(orders_r, pods_r, podswith_r, workstationassignment, ordersassignedto, orderassignmentbegan)

	model = Model(() -> Gurobi.Optimizer(GRB_ENV)) #Model(Gurobi.Optimizer)
	#set_optimizer_attribute(model, "TimeLimit", 600) # iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", 1) #ipoutputflag)

	#Variables
	@variable(model, h[m = orders_r, i = itemson[m], times], Bin)
	@variable(model, y[p = pods_r, a = podarcset[p]], Bin)
	@variable(model, v[orders_r, times], Bin)
	@variable(model, f[orders_r, times], Bin)
	@variable(model, g[orders_r, times], Bin)

	#Objective
	@objective(model, Max, sum(sum(h[m, i, last(times)] for i in itemson[m]) for m in orders_r) )

	#Linking constraints
	@constraint(model, podlink[m = orders_r, i = itemson[m], t = setdiff(times,times[1])], h[m,i,t] - h[m,i,t-tstep] <= sum(sum(sum(y[p,a] for a in intersect(podarcset[p], RA_plus[routenodes[workstationassignment[m],t2]])) for p in podswith_r[i,m]) for t2 in t:tstep_r:t+tstep-tstep_r))
	@constraint(model, podlinkzero[m = orders_r, i = itemson[m]], h[m,i,times[1]] <= sum(sum(sum(y[p,a] for a in intersect(podarcset[p], RA_plus[routenodes[workstationassignment[m],t2]])) for p in podswith_r[i,m]) for t2 in times[1]:tstep_r:times[1]+tstep-tstep_r) )
	@constraint(model, nondecreasing[m = orders_r, i = itemson[m], t = setdiff(times,times[1])], h[m,i,t] >= h[m,i,t-tstep] )

	#Pod movement constraints
	@constraint(model, podstartlocations[p = pods_r], sum(y[p,a] for a in intersect(podarcset[p], RA_plus[podstartnode_r[p]])) == 1)
	@constraint(model, podflowbalance[p = pods_r, n in setdiff(1:extendednumnodes_r, union(routenodes_end, podstartnode_r[p]))], sum(y[p,a] for a in intersect(podarcset[p], RA_minus[n])) - sum(y[p,a] for a in intersect(podarcset[p], RA_plus[n])) == 0)

	#Workstation capacity constraints
	@constraint(model, fconfig[m = orders_r, i = itemson[m], t = times], f[m,t] >= h[m,i,t])
	@constraint(model, gconfig[m = orders_r, i = itemson[m], t = times], g[m,t] <= h[m,i,t] )
	@constraint(model, vconfig[m = orders_r, t = times], v[m,t] == f[m,t] - g[m,t])
	@constraint(model, stationcapacity[w = workstations, t = times], sum(v[m,t] for m in ordersassignedto[w]) <= C[w] )
	@constraint(model, vearliest[m = orders_r, t = times], v[m,t] <= orderassignmentbegan[m,t])

	#Congestion
	@constraint(model, maxcongestion[l in intersections, t in -maxtraveltime:tstep_cong:horizon+maxtraveltime], sum(sum(y[p,a] for a in intersect(podarcset[p], congcontribarcs[l,t])) for p in pods_r) <= intersectionmaxpods[l])

	#====================================================#

	#Solve IP
	status = optimize!(model)

	if termination_status(model) == MOI.OPTIMAL
		solvetime = solve_time(model)
		obj = objective_value(model)
		#println("Total items picked = ", obj)
	elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
		#println("Time out in the subproblem")
		solvetime = solve_time(model)
		obj = objective_value(model)
		#println("Total items picked = ", obj)
	else
		println("Error in solving!")
		solvetime = solve_time(model)
		obj = 0
		#println("Total items picked = ", obj)
	end

	#====================================================#

	return obj, solvetime, h, y, v

end

#----------------------------------------------------------------------------------------------#

function solvedelayproblem(orders_r, pods_r, podswith_r, workstationassignment, ordersassignedto, orderassignmentbegan)

	model = Model(() -> Gurobi.Optimizer(GRB_ENV)) #Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", 60*60*12) # iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", 1) #ipoutputflag)

	#Variables
	@variable(model, y[p = pods_r, a = podarcset[p]], Bin)

	#Objective
	@objective(model, Min, sum(sum(sum(sum(poddelaypenalty[m,p,a] * y[p,a] for a in intersect(podarcset[p], RA_plus[routenodes[workstationassignment[m], t]])) for t in orderopentime[m]:tstep_r:horizon_r) for p in podsfor[m]) for m in orders_r) )

	#Pod delivery
	@constraint(model, vearliest[m = orders_r, p in podsfor[m]], sum(sum(y[p,a] for a in intersect(podarcset[p], A_space, RA_plus[routenodes[workstationassignment[m], t]])) for t in orderopentime[m]:tstep_r:horizon_r-tstep_r) >= 1)

	#Pod movement constraints
	@constraint(model, podstartlocations[p = pods_r], sum(y[p,a] for a in intersect(podarcset[p], RA_plus[podstartnode_r[p]])) == 1)
	@constraint(model, podflowbalance[p = pods_r, n in setdiff(1:extendednumnodes_r, union(routenodes_end, podstartnode_r[p]))], sum(y[p,a] for a in intersect(podarcset[p], RA_minus[n])) - sum(y[p,a] for a in intersect(podarcset[p], RA_plus[n])) == 0)

	#Congestion
	@constraint(model, maxcongestion[l in intersections, t in -maxtraveltime:tstep_cong:horizon_r+maxtraveltime], sum(sum(y[p,a] for a in intersect(podarcset[p], congcontribarcs[l,t])) for p in pods_r) <= intersectionmaxpods[l])

	#====================================================#

	#Solve IP
	#@objective(model, Min, sum(sum(sum(sum(poddelaypenalty[m,p,a] * y[p,a] for a in intersect(podarcset[p], RA_plus[routenodes[workstationassignment[m], t]])) for t in orderopentime[m]:tstep_r:horizon_r) for p in podsfor[m]) for m in orders_r) )
	status = optimize!(model)

	if termination_status(model) == MOI.OPTIMAL
		solvetime = solve_time(model)
		obj = objective_value(model)
		#println("Total items picked = ", obj)
	elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
		#println("Time out in the subproblem")
		solvetime = solve_time(model)
		obj = objective_value(model)
		#println("Total items picked = ", obj)
	else
		println("Error in solving!")
		solvetime = solve_time(model)
		obj = 0
		#println("Total items picked = ", obj)
	end

	#====================================================#

	return obj, solvetime, y

end

#=
totalinttimes, constrainedinttimes = 0, 0
totalslots, usedslots = 0, 0
for l in intersections, t in -maxtraveltime:tstep_cong:horizon+maxtraveltime
	total = 0
	for a in congcontribarcs[l,t], p in arcpodset[a]
		total += value(y[p,a]) 
	end
	totalinttimes += 1
	if total > intersectionmaxpods[l] - 1e-3
		constrainedinttimes += 1
	end
	totalslots += intersectionmaxpods[l]
	usedslots += total
end
println("Percent tight = ", constrainedinttimes / totalinttimes )
println("Percent used = ", usedslots / totalslots )

arcpodset = Dict()
for a in 1:numarcs_r
	arcpodset[a] = []
	for p in pods_r
		if a in podarcset[p]
			push!(arcpodset[a], p)
		end
	end
end

totaltraffic = Dict()
for l in intersections, t in -maxtraveltime:tstep_cong:horizon_r+maxtraveltime
	totaltraffic[l,t] = 0
end

for p in pods_r, a in podarcset[p]
	if value(y_nocong[p,a]) > 1e-4
		for (l,t) in passingintersectiontimes[a]
		totaltraffic[l,t] += value(y_nocong[p,a]) 
		end
	end
end

for l in intersections, t in -maxtraveltime:tstep_cong:horizon+maxtraveltime
	if totaltraffic[l,t] > intersectionmaxpods[l]
	println("totaltraffic[$l,$t] = ", totaltraffic[l,t], " <= ", intersectionmaxpods[l])
	end
end

	
for m in orders_r, p in podsfor[m], t in orderopentime[m]:tstep_r:horizon_r
	for a in intersect(podarcset[p], RA_plus[routenodes[workstationassignment[m], t]]) 
		if (value(y[p,a]) > 1e-4) & ( poddelaypenalty[m,p,a] > 1e-4)
			println("Pod $p for order $m delivered at time $t (planned = ", orderopentime[m], " - ", orderclosetime[m], ") incurring delay = ", poddelaypenalty[m,p,a])
		end
	end
end


#First solution --> delay = 75 (dynamic lol)
Pod 51 for order 319 delivered at time 180.0 (planned = 120 - 150) incurring delay = 1.0
Pod 334 for order 319 delivered at time 750.0 (planned = 120 - 150) incurring delay = 20.0
Pod 545 for order 319 delivered at time 420.0 (planned = 120 - 150) incurring delay = 9.0
Pod 200 for order 497 delivered at time 660.0 (planned = 420 - 600) incurring delay = 2.0
Pod 385 for order 1496 delivered at time 750.0 (planned = 510 - 540) incurring delay = 7.0
Pod 74 for order 1023 delivered at time 660.0 (planned = 540 - 570) incurring delay = 3.0
Pod 113 for order 652 delivered at time 840.0 (planned = 690 - 720) incurring delay = 4.0
Pod 411 for order 652 delivered at time 750.0 (planned = 690 - 720) incurring delay = 1.0
Pod 471 for order 469 delivered at time 660.0 (planned = 0 - 150) incurring delay = 17.0
Pod 132 for order 42 delivered at time 600.0 (planned = 330 - 570) incurring delay = 1.0
Pod 334 for order 42 delivered at time 690.0 (planned = 330 - 570) incurring delay = 4.0
Pod 321 for order 391 delivered at time 390.0 (planned = 0 - 300) incurring delay = 3.0
Pod 653 for order 100 delivered at time 600.0 (planned = 300 - 510) incurring delay = 3.0

#Other solution --> delay = 61 (greedy)
Pod 269 for order 637 delivered at time 690.0 (planned = 0 - 270) incurring delay = 14.0
Pod 28 for order 469 delivered at time 840.0 (planned = 330 - 630) incurring delay = 7.0
Pod 545 for order 244 delivered at time 630.0 (planned = 390 - 600) incurring delay = 1.0
Pod 583 for order 889 delivered at time 660.0 (planned = 60 - 180) incurring delay = 16.0
Pod 647 for order 1704 delivered at time 330.0 (planned = 60 - 150) incurring delay = 6.0
Pod 150 for order 259 delivered at time 540.0 (planned = 120 - 480) incurring delay = 2.0
Pod 172 for order 554 delivered at time 660.0 (planned = 150 - 240) incurring delay = 14.0
Pod 465 for order 1777 delivered at time 660.0 (planned = 540 - 630) incurring delay = 1.0


usedtime = 0
delaytime = 0
itemtime, podtime = 0, 0 
for m in orders_r, p in podsfor[m], t in orderopentime[m]:tstep_r:horizon
	for a in intersect(podarcset[p], RA_plus[routenodes[workstationassignment[m], t]]) 
		if (value(y[p,a]) > 1e-4) && (a in A_space)
			usedtime += itemprocesstime * length(itemsfrom[p,m]) + podprocesstime
			itemtime += itemprocesstime * length(itemsfrom[p,m])
			podtime += podprocesstime
			delaytime += (t - orderopentime[m])
			if (t - orderopentime[m]) > 1e-4
				println("Pod $p for order $m ==> ", t, " - ", orderopentime[m])
			end
		end
	end
end
idletime = horizon*numworkstations - usedtime
itemtime, podtime, delaytime
println("Delay per pod = ", delaytime / (podtime/podprocesstime))



for m in orders_r
	println("$m --> ", orderopentime[m], " - ", orderclosetime[m])
	for w in workstations, t in 0:tstep:horizon_r, (m2,i,p) in newpicklist[w,t]
		if m == m2
			println("     ($m2, $i, $p, $w, $t)")
		end
	end
end

for l in queueintersections, t in -maxtraveltime:tstep_cong:horizon_r+maxtraveltime
	traffic = 0
	for p in pods_r, a in intersect(podarcset[p], congcontribarcs[l,t])
		traffic += value(y_r[p,a])
	end
	if traffic >= intersectionmaxpods[l] -100
		println(traffic, " <= ", intersectionmaxpods[l])
	end
end


=#

#------------------------------------------------------------------------------------------------------#

function writeroutingmetrics(filename, itempodpicklist, y)

	#----------------------------------------- PLAN -----------------------------------------#

	#Throughput, idle time and utilization
	usedtime, itemtime, podtime, throughput_plan = 0, 0, 0, 0
	for w in workstations, t in 0:tstep:horizon
		usedpods = []
		for (m,i,p) in itempodpicklist[w,t]
			usedpods = union(usedpods, p)
			usedtime += itemprocesstime
			itemtime += itemprocesstime
			throughput_plan += 1
		end
		for p in usedpods
			usedtime += podprocesstime
			podtime += podprocesstime
		end
	end
	idle_plan = (horizon+tstep)*numworkstations - usedtime
	utilization_plan = 1 - (idle_plan / ((horizon+tstep)*numworkstations))
	println("Thoughput = ", throughput_plan)
	println("Idle time = ", idle_plan)
	println("Utilizatn = ", utilization_plan)

	#Find original delivery times
	planneddeliverytime = Dict()
	for w in workstations, t in 0:tstep_r:horizon, (m,i,p) in itempodpicklist[w,t]
		planneddeliverytime[m,i,p] = t
	end
	plannedclosetime = Dict()
	for m in orders_r
		plannedclosetime[m] = -1
	end
	for (m,i,p) in keys(planneddeliverytime)
		plannedclosetime[m] = max(plannedclosetime[m], planneddeliverytime[m,i,p])
	end

	#---------------------------------------- ACTUAL ----------------------------------------#

	#Get an updated itempodpicklist
	newpicklist = Dict()
	for w in workstations, t in 0:tstep_r:horizon_r
		newpicklist[w,t] = []
		for m in ordersassignedto[w], p in podsfor[m]
			if t >= orderopentime[m]
				for a in intersect(podarcset[p], RA_plus[routenodes[workstationassignment[m], t]]) 
					if (value(y[p,a]) > 1e-4) && (a in A_space)
						for i in itemsfrom[p,m]
							push!(newpicklist[w,t], (m,i,p))
						end
					end
				end
			end
		end 
		newpicklist[w,t] = unique(newpicklist[w,t])
	end
	
	#Throughput, idle time, utilization, and delay
	usedtime, numitems, numpods, delaytime_orig, delaytime, delaytime_order, throughput_actual = 0, 0, 0, 0, 0, 0, 0
	for w in workstations, t in 0:tstep:horizon
		usedpods = []
		for (m,i,p) in newpicklist[w,t]
			#push!(usedpods, (p, max(0,t - planneddeliverytime[m,i,p])))
			push!(usedpods, (p, max(0,t - planneddeliverytime[m,i,p]), t - planneddeliverytime[m,i,p], max(0,t - plannedclosetime[m])))
			usedtime += itemprocesstime
			throughput_actual += 1
			numitems += 1
		end
		usedpods = unique(usedpods)
		for p in usedpods
			usedtime += podprocesstime
			numpods += 1
			delaytime_orig += p[2]
			delaytime += p[3]
			delaytime_order += p[4]
			if p[2] > 1e-4
			end
		end
	end
	idle_actual = (horizon+tstep)*numworkstations - usedtime
	utilization_actual = 1 - idle_actual / ((horizon+tstep)*numworkstations)
	total_delay_actual = delaytime 
	pod_delay_actual = delaytime / numpods
	total_order_delay_actual = delaytime_order 
	order_pod_delay_actual = delaytime_order / numpods
	orig_pod_delay_actual = delaytime_orig / numpods
	println("Thoughput = ", throughput_actual)
	println("Idle time = ", idle_actual)
	println("Utilizatn = ", utilization_actual)	
	println("Ttl delay = ", total_delay_actual)	
	println("Pod delay (negatives) = ", pod_delay_actual)	
	println("Pod delay (order) = ", order_pod_delay_actual)	
	println("Pod delay (orig) = ", orig_pod_delay_actual)	

	#---------------------------------------- OUTPUT ----------------------------------------#

	df = DataFrame(
		run_id = [run_id],
		instance_id = [instance_id],
		warehouse_id = [warehouse_id],
		method = [methodname],
		horizon = [horizon_r-tstep_r],
		throughput_plan = [throughput_plan],
		idle_plan = [idle_plan],
		utilization_plan = [utilization_plan],
		throughput_actual = [throughput_actual],
		idle_actual = [idle_actual],
		utilization_actual = [utilization_actual],
		total_delay_actual = [total_delay_actual],
		pod_delay_actual = [pod_delay_actual]
		)

	CSV.write(filename, df)

end