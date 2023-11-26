
function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end

#-----------------------------------------------------------------------------------#

function initializeblanksolution()

	podsworked, ordersworking, itemspicked = Dict(), Dict() , Dict()
	for w in workstations, t in times
		podsworked[w,t] = []
		ordersworking[w,t] = []
		itemspicked[w,t] = []
	end

	remaininginventory = deepcopy(inventory)
	podswith_greedy = deepcopy(podswith)

	podawayfromhome = Dict()
	podjourneys = Dict()
	itemsremoved = Dict()

	for p in allpods
		podjourneys[p] = []
		podawayfromhome[p] = []
		for t in times
			itemsremoved[p,t] = []
		end
	end

	congestionat = Dict()
	for i in reducedintersections, t in 0:congestiontstep:horizon
		congestionat[i,t] = 0
	end

	return podsworked, ordersworking, remaininginventory, podswith_greedy, podawayfromhome, podjourneys, itemspicked, itemsremoved, congestionat

end

#-----------------------------------------------------------------------------------#

function prioritizelargerorders(orderlist)

	#sortedpairings = sort!([pair for pair in itemson], by = x -> length(x[2]), rev=true)
	#sortedorders = [item[1] for item in sortedpairings]

	sortedorders = reverse(sort(orderlist, by=x->length(itemson[x])))

	return sortedorders

end

#-----------------------------------------------------------------------------------#

function prioritizeorders(orderlist, orderprioritization)

	if orderprioritization == 0
		sortedorders = shuffle(orderlist)
	elseif orderprioritization == 1
		sortedorders =  prioritizelargerorders(orderlist)
	elseif orderprioritization == -1
		sortedorders = shuffle(orderlist)
	else
		println("Order prioritization flag not recognized")
	end

	return sortedorders

end

#-----------------------------------------------------------------------------------#

function checkpodoverlap(m, podswith_greedy, podsworked)

	relevantpods = []
	for i in itemson[m]
		relevantpods = union(relevantpods, podswith_greedy[i])
	end

	overlapopportunities = Dict()
	for w in workstations, t in times
		overlapopportunities[w,t] = intersect(relevantpods, podsworked[w,t])
	end

	workstationoverlap = Dict()
	for w in workstations
		workstationoverlap[w] = sum(length(overlapopportunities[w,t]) for t in times)
	end

	return overlapopportunities, workstationoverlap, relevantpods

end

#-----------------------------------------------------------------------------------#

function manhattandist(l,w)

	dist = abs(loccoords[l,1] - loccoords[w,1]) + abs(loccoords[l,2] - loccoords[w,2])

	return dist

end

#-----------------------------------------------------------------------------------#

function orderclosestworkstations(m, podswith_greedy)

	totalworkstationdistance = Dict()

	for w in workstations
		totalworkstationdistance[w] = 0
		for i in itemson[m]
			closestpod = 10000000
			for p in podswith_greedy[i]
				closestpod = min(closestpod, manhattandist(podstorageloc[p],w))
			end
			totalworkstationdistance[w] += closestpod
		end
	end

	return totalworkstationdistance

end

#-----------------------------------------------------------------------------------#

#capacitybuffer=capacitybuffer2
function fitorderin(m, w, overlapopportunities, capacitybuffer, ordersworking, podsworked, itemspicked)

	model = Model(() -> Gurobi.Optimizer(GRB_ENV)) #Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "OutputFlag", 0)
	set_optimizer_attribute(model, "TimeLimit", 10)

	@variable(model, pick[times] >= 0, Int)
	@variable(model, firsti[times], Bin)
	@variable(model, lasti[times], Bin)

	@objective(model, Min, sum((firsti[t] - lasti[t]) for t in times))

	@constraint(model, sum(pick[t] for t in times) == length(itemson[m]))
	@constraint(model, [t in times], firsti[t] >= sum(pick[t2] for t2 in 0:tstep:t)/ length(itemson[m]))
	@constraint(model, [t in times], lasti[t] <= sum(pick[t2] for t2 in 0:tstep:t)/ length(itemson[m]))
	@constraint(model, [t in times], (itemprocesstime + podprocesstime - podprocesstime * min(1,length(overlapopportunities[w,t]))) * pick[t] <= max(0, tstep - podprocesstime * length(podsworked[w,t]) - itemprocesstime * length(itemspicked[w,t]) - capacitybuffer) )
	if length(itemson[m]) >= 2
		@constraint(model, [t in times], firsti[t] - lasti[t] + length(ordersworking[w,t]) <= workstationordercapacity)
	end

	optimize!(model)

	if (termination_status(model) == MOI.OPTIMAL) || ((termination_status(model) == MOI.TIME_LIMIT) & (primal_status(model) == MOI.FEASIBLE_POINT))
		orderstarttime, orderendtime = 0, horizon
		feasible_flag = 1
		loop2_flag = 0
		for t in times
			if (loop2_flag == 0) & (getvalue(firsti[t]) - getvalue(lasti[t]) > 0.01)
				orderstarttime = t
				loop2_flag = 1
			elseif (loop2_flag == 1) & (getvalue(lasti[t]) > 0.01)
				orderendtime = t
				break
			end
		end
		if orderstarttime == 0
			orderstarttime = dummystarttime
		end
		if orderendtime == horizon
			orderendtime = dummyendtime
		end
	else
		feasible_flag = 0
		orderstarttime, orderendtime = -1, -1
	end

	return feasible_flag, orderstarttime, orderendtime

end

#-----------------------------------------------------------------------------------#

#Tries to make official assignment for order m at workstation w
#m, w, orderstarttime, orderendtime, sp_relevantpods, sp_relevantpods_t, podswith_greedy, podsworked, itemspicked, congestionat, itemsremoved, remaininginventory, ordersworking, podstarttime, podendtime, podawayfromhome, stationassign, ordersassignedto = m, w, orderstarttime, orderendtime, sp_relevantpods, sp_relevantpods_t, podswith_greedy, podsworked, itemspicked, congestionat, itemsremoved, remaininginventory, ordersworking, podstarttime, podendtime, podawayfromhome, stationassign, ordersassignedto

function findfeasibleorderassignment(m, w, orderstarttime, orderendtime, sp_relevantpods, sp_relevantpods_t, podswith_greedy, podsworked, itemspicked, congestionat, itemsremoved, remaininginventory, ordersworking, podstarttime, podendtime, podawayfromhome, stationassign, ordersassignedto)
	
	sp_times = orderstarttime:tstep:orderendtime
	sp_times_reg = max(0,orderstarttime):tstep:min(horizon,orderendtime)

	model = Model(() -> Gurobi.Optimizer(GRB_ENV)) #Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", 10) # iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", 0)

	#Variables
	@variable(model, h_m[i = itemson[m], p = intersect(sp_relevantpods, podswith_greedy[i]), podstarttime[p]:tstep:podendtime[p]], Bin)
	@variable(model, y_m[p = sp_relevantpods, t = podstarttime[p]:tstep:podendtime[p]], Bin) #if 1 ==> t is arrival time of pod 
	@variable(model, z_m[p = sp_relevantpods, t = podstarttime[p]:tstep:podendtime[p]], Bin) #if 1 ==> t is departure time of pod 
	@variable(model, w_m[p = sp_relevantpods, t = podstarttime[p]:tstep:podendtime[p]], Bin) #if 1 ==> pod sits in queue at time t
	@variable(model, m_incomplete, Bin)

	#Objective
	@objective(model, Max, sum(sum(sum(h_m[i, p, t] for t in podstarttime[p]:tstep:podendtime[p]) for p in intersect(sp_relevantpods, podswith_greedy[i])) for i in itemson[m]) - 0.1 * m_incomplete - 0.001 * sum(sum(y_m[p,t] + z_m[p,t] + w_m[p,t] for t in podstarttime[p]:tstep:podendtime[p]) for p in sp_relevantpods) )

	#Order delivery
	if orderendtime >= horizon
		@constraint(model, orderdelivery_max[i = itemson[m]], sum(sum(h_m[i, p, t] for t in podstarttime[p]:tstep:podendtime[p]) for p in intersect(sp_relevantpods, podswith_greedy[i])) <= 1)
		@constraint(model, orderdelivery_min[i = itemson[m]], sum(sum(h_m[i, p, t] for t in podstarttime[p]:tstep:podendtime[p]) for p in intersect(sp_relevantpods, podswith_greedy[i])) >= 1 - m_incomplete)
	else
		@constraint(model, orderdelivery_max[i = itemson[m]], sum(sum(h_m[i, p, t] for t in podstarttime[p]:tstep:podendtime[p]) for p in intersect(sp_relevantpods, podswith_greedy[i])) == 1)
	end

	#Item-pod link
	#@constraint(model, podonlyifitem[i = itemson[m], p = intersect(sp_relevantpods, podswith_greedy[i]), t = podstarttime[p]:tstep:podendtime[p], t2 = t:tstep:podendtime[p]], y_m[p,t] <= h_m[i,p,t2])
	@constraint(model, podarrives[i = itemson[m], p = intersect(sp_relevantpods, podswith_greedy[i]), t = podstarttime[p]:tstep:podendtime[p]], sum(y_m[p,t2] for t2 in podstarttime[p]:tstep:t) >= h_m[i,p,t])
	@constraint(model, poddeparts[i = itemson[m], p = intersect(sp_relevantpods, podswith_greedy[i]), t = podstarttime[p]:tstep:podendtime[p]], z_m[p,t] >= h_m[i,p,t])
	@constraint(model, podwaits[p = sp_relevantpods, t = podstarttime[p]:tstep:podendtime[p]], w_m[p,t] >= sum(y_m[p,t2] for t2 in podstarttime[p]:tstep:t) - sum(z_m[p,t2] for t2 in podstarttime[p]:tstep:t))

	#Workstation throughput constraints
	@constraint(model, stationthroughput[t = sp_times_reg], itemprocesstime * sum(sum(h_m[i,p,t] for p in intersect(sp_relevantpods, podswith_greedy[i]) if podstarttime[p] <= t <= podendtime[p]) for i in itemson[m]) + podprocesstime * sum(z_m[p,t] for p in sp_relevantpods_t[t]) <= tstep - itemprocesstime * length(itemspicked[w,t]) - podprocesstime * length(podsworked[w,t]))

	#Congestion
	#=
	for i in setdiff(reducedintersections,queueintersections), t in max(0,orderstarttime):congestiontstep:min(horizon-tstep,orderendtime)
		@constraint(model, 
			#sum(sum(congestioncontribution[arcs[nodes[w, t2], nodes[podstorageloc[p], t2 + arclength[w,podstorageloc[p]]]], i, t] * z_m[p,t2] for t2 in podstarttime[p]:tstep:podendtime[p] if arcs[nodes[w, t2], nodes[podstorageloc[p], t2 + arclength[w,podstorageloc[p]]]] in congestionarcs[i,t]) for p in sp_relevantpods) +
			#sum(sum(congestioncontribution[arcs[nodes[podstorageloc[p], t2 - arclength[podstorageloc[p],w]], nodes[w, t2]], i, t] * y_m[p,t2] for t2 in podstarttime[p]:tstep:podendtime[p] if arcs[nodes[podstorageloc[p], t2 - arclength[podstorageloc[p],w]], nodes[w, t2]] in congestionarcs[i,t]) for p in sp_relevantpods)
			sum(sum(sum(congestioncontribution[a, i, t] * z_m[p,t2] for a in intersect(congestionarcs[i,t], sp_zarc[w,p,t2])) for t2 in podstarttime[p]:tstep:podendtime[p]) for p in sp_relevantpods)
			+ sum(sum(sum(congestioncontribution[a, i, t] * y_m[p,t2] for a in intersect(congestionarcs[i,t], sp_yarc[w,p,t2])) for t2 in podstarttime[p]:tstep:podendtime[p]) for p in sp_relevantpods)
			 <= max(0, intersectionmaxpods[i] - congestionat[i,t]))
	end
	
	#Workstation queue congestion
	i_q = maploctointersection[w]
	for t in max(0,orderstarttime):congestiontstep:min(horizon-tstep,orderendtime)
		t_q = floor(t/tstep) * tstep
		if t_q <= horizon - tstep
			a_q = arcs[nodes[w,t_q], nodes[w,t_q+tstep]]
			@constraint(model, 
				sum(sum(sum(congestioncontribution[a, i_q, t] * z_m[p,t2] for a in intersect(congestionarcs[i_q,t], sp_zarc[w,p,t2])) for t2 in podstarttime[p]:tstep:podendtime[p]) for p in sp_relevantpods)
				+ sum(sum(sum(congestioncontribution[a, i_q, t] * y_m[p,t2] for a in intersect(congestionarcs[i_q,t], sp_yarc[w,p,t2])) for t2 in podstarttime[p]:tstep:podendtime[p]) for p in sp_relevantpods)
				+ sum(congestioncontribution[a_q, i_q, t] * w_m[p, t_q] for p in sp_relevantpods if t_q in podstarttime[p]:tstep:podendtime[p])
				<= max(0, intersectionmaxpods[i_q] - congestionat[i_q,t]))
		end
	end
	=#

	#====================================================#

	#Solve IP
	status = optimize!(model)

	if termination_status(model) == MOI.OPTIMAL
		solvetime = solve_time(model)
		obj = 0
		for i in itemson[m], p in intersect(sp_relevantpods, podswith_greedy[i]), t in podstarttime[p]:tstep:podendtime[p]
			obj += getvalue(h_m[i, p, t])
		end
		if obj > 0
			println("Found assignment")
			feasibleflag = 1
		else
			println("Order not assigned")
			feasibleflag = 0
		end
	elseif (termination_status(model) == MOI.TIME_LIMIT) & (primal_status(model) == MOI.FEASIBLE_POINT)
		solvetime = solve_time(model)
		obj = 0
		for i in itemson[m], p in intersect(sp_relevantpods, podswith_greedy[i]), t in podstarttime[p]:tstep:podendtime[p]
			obj += getvalue(h_m[i, p, t])
		end
		if obj > 0
			println("Found assignment (timeout)")
			feasibleflag = 1
		else
			println("Order not assigned")
			feasibleflag = 0
		end
	else
		println("Order not assigned")
		solvetime = solve_time(model)
		obj = 0
		feasibleflag = 0
	end

	#====================================================#

	#If we found a feasible assignment, update the solution
	if feasibleflag == 1

		earliestt, latestt = min(horizon, orderendtime), max(0, orderstarttime)
		for i in itemson[m]
			for p in intersect(sp_relevantpods, podswith_greedy[i]), t in podstarttime[p]:tstep:podendtime[p]
				if getvalue(h_m[i,p,t]) > 0.01
					push!(itemsremoved[p, t], i)
					push!(itemspicked[w, t], (m,i,p))
					remaininginventory[i, p] -= 1
					earliestt = min(earliestt, t)
					latestt = max(latestt, t) #Here is the fix for the greedy algorithm!
					podsworked[w,t] = union(podsworked[w,t], p)
					remove!(podswith_greedy[i], p)
					push!(itempodpicklist[w, t], (m, i, p))
					for t_prime in t:tstep:horizon
						h_currsol[m, i, p, w, t_prime] = 1
					end
					break
				end
			end
		end

		for t in earliestt:tstep:latestt
			push!(ordersworking[w, t], m)
			v_currsol[m, w, t] = 1
		end

		for p in sp_relevantpods
			#podtostationtime = tstep * ceil(traveltimeraw[maploctointersection[podstorageloc[p]], maploctointersection[w]] / tstep)
			podtostationtime = arclength[podstorageloc[p], w]
			awaystart, awayend = 0, 0
			podarrives = 0
			for t in podstarttime[p]:tstep:podendtime[p]
				if getvalue(y_m[p,t]) > 0.01
					awaystart += max(0, t - podtostationtime)
					podarrives += t
				end
				if getvalue(z_m[p,t]) > 0.01
					awayend += min(horizon, t + podtostationtime)
					for t_p in awaystart:tstep:awayend
						push!(podawayfromhome[p], t_p)
					end

					#Pod trip to the station
					if t - podtostationtime < 0
						y_currsol[p, extendedarcs[extendednodes[podstorageloc[p], dummystarttime], extendednodes[w, podarrives]]] = 1
					else
						y_currsol[p, arcs[nodes[podstorageloc[p], awaystart], nodes[w, awaystart + podtostationtime]]] = 1
					end		
										
					#Pod trip back to storage region
					if t + podtostationtime > horizon
						y_currsol[p, extendedarcs[extendednodes[w, t], extendednodes[podstorageloc[p], dummyendtime]]] = 1	
					else 
						y_currsol[p, arcs[nodes[w, awayend - podtostationtime], nodes[podstorageloc[p], awayend]]] = 1	
					end	

					#Pod stay in the workstation queue
					for t_p in awaystart+podtostationtime:tstep:min(horizon-tstep, t-tstep)
						y_currsol[p,arcs[nodes[w, t_p], nodes[w, t_p + tstep]]] = 1
					end

					#Pod no longer stationed at home
					for t_p in awaystart:tstep:min(dummyendtime-tstep, awayend - tstep)
						y_currsol[p, arcs[nodes[podstorageloc[p], t_p], nodes[podstorageloc[p], t_p + tstep]]] = 0
					end
					break
				end
			end
		end

		#AUDIT THE CONGESTIONAT UPDATE!

		for i in setdiff(reducedintersections, queueintersections), t in max(0, orderstarttime):congestiontstep:min(horizon, orderendtime)
			for (s,t2) in congestionpairsy[w,i,t] 
				if (t2 in sp_times_reg)
					for p in intersect(sp_relevantpods_t[t2], podsstoredat[s])
						congestionat[i,t] += congestioncontributiony[podstorageloc[p],w,i,t2,t] * getvalue(y_m[p,t2])
					end
				end
			end
			for (s,t2) in congestionpairsz[w,i,t]
				if (t2 in sp_times_reg)
					for p in intersect(sp_relevantpods_t[t2], podsstoredat[s])
						congestionat[i,t] += congestioncontributionz[podstorageloc[p],w,i,t2,t] * getvalue(z_m[p,t2])
					end
				end
			end
			#congestionat[i,t] += sum(sum(congestioncontributiony[podstorageloc[p],w,i,t2,t] * getvalue(y_m[p,t2]) for p in sp_relevantpods_t[t2]) for t2 in sp_times if sp_relevantpods_t[t2] != [])
			#congestionat[i,t] += sum(sum(congestioncontributionz[podstorageloc[p],w,i,t2,t] * getvalue(z_m[p,t2]) for p in sp_relevantpods_t[t2]) for t2 in sp_times if sp_relevantpods_t[t2] != [])
		end

		for t in max(0, orderstarttime):congestiontstep:min(horizon, orderendtime)
			i = maploctointersection[w]
			for (s,t2) in congestionpairsy[w,i,t] 
				if (t2 in sp_times_reg)
					for p in intersect(sp_relevantpods_t[t2], podsstoredat[s])
						congestionat[i,t] += congestioncontributiony[podstorageloc[p],w,i,t2,t] * getvalue(y_m[p,t2])
					end
				end
			end
			for (s,t2) in congestionpairsz[w,i,t]
				if (t2 in sp_times_reg)
					for p in intersect(sp_relevantpods_t[t2], podsstoredat[s])
						congestionat[i,t] += congestioncontributionz[podstorageloc[p],w,i,t2,t] * getvalue(z_m[p,t2])
					end
				end
			end
			t_q = floor(t/tstep) * tstep
			if t_q <= horizon - tstep
				a_q = arcs[nodes[w,t_q], nodes[w,t_q+tstep]]
			 	for p in sp_relevantpods_t[t_q]
			 		congestionat[i,t] += congestioncontribution[a_q, i, t] * getvalue(w_m[p,t_q])
			 	end
			end
		end

		remove!(unassignedorders, m)
		stationassign[m] = w
		push!(ordersassignedto[w], m)

		#println("    Instance updated")
	end

	#====================================================#

	return obj, solvetime, h_m, y_m, z_m, feasibleflag, podswith_greedy, podsworked, itemspicked, congestionat, itemsremoved, remaininginventory, ordersworking, podawayfromhome, stationassign, ordersassignedto

end

#-----------------------------------------------------------------------------------#

function getpodstarttimes(relevantpods, orderstarttime, orderendtime, podawayfromhome, w)

	podstarttime, podendtime = Dict(), Dict()
	sp_relevantpods = []
	sp_relevantpods_t = Dict()
	for t in orderstarttime:congestiontstep:orderendtime
		sp_relevantpods_t[t] = []
	end

	for p in relevantpods
		podtostationtime = arclength[podstorageloc[p], w]
		podfreetimes = setdiff(orderstarttime:tstep:orderendtime, podawayfromhome[p])
		if (podfreetimes != []) && (length(podfreetimes) == (last(podfreetimes) - podfreetimes[1]) / tstep + 1)
			push!(sp_relevantpods, p)
			if podfreetimes[1] == dummystarttime
				podstarttime[p] = 0
			else
				podstarttime[p] = max(podfreetimes[1] + podtostationtime, orderstarttime)
			end
			#### CHANGE ####
			if last(podfreetimes) == dummyendtime
				podendtime[p] = horizon
			else
				podendtime[p] = min(last(podfreetimes) - podtostationtime + tstep, orderendtime)
			end
			for t2 in podstarttime[p]:congestiontstep:podendtime[p]
				push!(sp_relevantpods_t[t2], p)
			end
		end
	end

	return sp_relevantpods, sp_relevantpods_t, podstarttime, podendtime

end

#-----------------------------------------------------------------------------------#

#Input: instance and parameters
#Output: a greedy heuristic solution of the instance
#Proceed:
function getgreedysolution(orders, orderprioritization)

	stationassign, ordersassignedto = Dict(), Dict()
	for w in workstations
		ordersassignedto[w] = []
	end

	podsworked, ordersworking, remaininginventory, podswith_greedy, podawayfromhome, podjourneys, itemspicked, itemsremoved, congestionat = initializeblanksolution()

	sortedorders = prioritizeorders(orders, orderprioritization)

	totalobj = 0

	for m in sortedorders
		println("----- ORDER $m (", length(itemson[m]), " items) -----")
		
		overlapopportunities, workstationoverlap, relevantpods = checkpodoverlap(m, podswith_greedy, podsworked)
		totalworkstationdistance = orderclosestworkstations(m, podswith_greedy)
		sortedworkstationtups = sort!([(w, workstationoverlap[w], totalworkstationdistance[w]) for w in workstations], by = x -> 1000*x[2] - 0.001*x[3], rev=true)
		sortedworkstations = [item[1] for item in sortedworkstationtups]
		
		for w in sortedworkstations
			capacitybuffer2 = 5*min(floor(length(itemson[m])/2),3)
			feasible_flag, orderstarttime, orderendtime = fitorderin(m, w, overlapopportunities, capacitybuffer2, ordersworking, podsworked, itemspicked)
			if feasible_flag == 1
				sp_times = [t2 for t2 in orderstarttime:tstep:orderendtime]
				sp_relevantpods, sp_relevantpods_t, podstarttime, podendtime = getpodstarttimes(relevantpods, orderstarttime, orderendtime, podawayfromhome, w)
				assign_obj, assign_solvetime, assign_h, assign_y, assign_z, assignmentsuccessful, podswith_greedy, podsworked, itemspicked, congestionat, itemsremoved, remaininginventory, ordersworking, podawayfromhome, stationassign, ordersassignedto = findfeasibleorderassignment(m, w, orderstarttime, orderendtime, sp_relevantpods, sp_relevantpods_t, podswith_greedy, podsworked, itemspicked, congestionat, itemsremoved, remaininginventory, ordersworking, podstarttime, podendtime, podawayfromhome, stationassign, ordersassignedto)
				if assignmentsuccessful == 1
					totalobj += floor(assign_obj)
					break
				end
			end
		end
	end
				
	totaltime = horizon*numworkstations
	usedtime = 0
	for t in times, w in workstations
		usedtime += itemprocesstime * length(itemspicked[w,t]) + podprocesstime * length(podsworked[w,t])
	end
	println("Final utilization = ", round(usedtime / totaltime * 100, digits = 2), "%")

	y_currpath = Dict()
	for p in pods
		y_currpath[p] = []
		for a in podarcset[p]
			if y_currsol[p,a] == 1
				push!(y_currpath[p], a)
			end
		end
	end

	println("Passed")

	return stationassign, ordersassignedto, y_currpath, ordersworking, totalobj, itemspicked, usedtime

	#return podsworked, ordersopen, remaininginventory, podswith_greedy, podawayfromhome, podjourneys, itemspicked, itemsremoved, congestionat

end

#-----------------------------------------------------------------------------------#

function printcurrentsolution()

	for w in workstations
		println("----- STATION $w -----")
		for t in times
			println(" t = $t --> ", length(itemspicked[w,t]), " items picked from ", length(podsworked[w,t]), " pods")
		end
	end

end
