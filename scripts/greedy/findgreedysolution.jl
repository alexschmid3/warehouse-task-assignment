
include("updategreedysolution.jl")

#-----------------------------------------------------------------------------------#

function initializeblanksolution(podswith, allpods)

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

    currsol_greedysets = (remaininginventory=remaininginventory, podswith_greedy=podswith_greedy, podawayfromhome=podawayfromhome, podjourneys=podjourneys, itemsremoved=itemsremoved)

	return currsol_greedysets

end

#-----------------------------------------------------------------------------------#

function prioritizelargerorders(orderlist)

	sortedorders = reverse(sort(orderlist, by=x->length(itemson[x])))

	return sortedorders

end

#-----------------------------------------------------------------------------------#

function checkpodoverlap(m, workstations, currsol_greedysets, currsol)

	relevantpods = []
	for i in itemson[m]
		relevantpods = union(relevantpods, currsol_greedysets.podswith_greedy[i])
	end

	overlapopportunities = Dict()
	for w in workstations, t in times
		overlapopportunities[w,t] = intersect(relevantpods, currsol.podsworkedat[w,t])
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

function orderclosestworkstations(m, workstations, currsol_greedysets)

	totalworkstationdistance = Dict()

	for w in workstations
		totalworkstationdistance[w] = 0
		for i in itemson[m]
			closestpod = 10000000
			for p in currsol_greedysets.podswith_greedy[i]
				closestpod = min(closestpod, manhattandist(podstorageloc[p],w))
			end
			totalworkstationdistance[w] += closestpod
		end
	end

	return totalworkstationdistance

end

#-----------------------------------------------------------------------------------#

function fitorderin(m, w, overlapopportunities, capacitybuffer, currsol)

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
	@constraint(model, [t in times], (itemprocesstime + podprocesstime - podprocesstime * min(1,length(overlapopportunities[w,t]))) * pick[t] <= max(0, tstep - podprocesstime * length(currsol.podsworkedat[w,t]) - itemprocesstime * length(currsol.itempodpicklist[w,t]) - capacitybuffer) )
	if length(itemson[m]) >= 2
		@constraint(model, [t in times], firsti[t] - lasti[t] + length(currsol.ordersopen[w,t]) <= workstationordercapacity)
	end

	optimize!(model)

	if (termination_status(model) == MOI.OPTIMAL) || ((termination_status(model) == MOI.TIME_LIMIT) & (primal_status(model) == MOI.FEASIBLE_POINT))
		orderstarttime, orderendtime = 0, horizon
		feasible_flag = 1
		loop2_flag = 0
		for t in times
			if (loop2_flag == 0) & (value(firsti[t]) - value(lasti[t]) > 0.01)
				orderstarttime = t
				loop2_flag = 1
			elseif (loop2_flag == 1) & (value(lasti[t]) > 0.01)
				orderendtime = t
				break
			end
		end
		#if orderstarttime == 0
		#	orderstarttime = dummystarttime
		#end
		#if orderendtime == horizon
		#	orderendtime = dummyendtime
		#end
	else
		feasible_flag = 0
		orderstarttime, orderendtime = -1, -1
	end

	return feasible_flag, orderstarttime, orderendtime

end

#-----------------------------------------------------------------------------------#

#Tries to make official assignment for order m at workstation w
function findfeasibleorderassignment(gp, currpartition, currsol, currsol_greedysets)
	
    sp_times = gp.orderstarttime:tstep:gp.orderendtime
	sp_times_reg = max(0,gp.orderstarttime):tstep:min(horizon,gp.orderendtime)

	model = Model(() -> Gurobi.Optimizer(GRB_ENV)) #Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", 10) # iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", 0)

	#Variables
	@variable(model, h_m[i = itemson[gp.m], p = intersect(gp.relevantpods, currsol_greedysets.podswith_greedy[i]), gp.podstarttime[p]:tstep:gp.podendtime[p]], Bin)
	@variable(model, y_m[p = gp.relevantpods, t = gp.podstarttime[p]:tstep:gp.podendtime[p]], Bin) #if 1 ==> t is arrival time of pod 
	@variable(model, z_m[p = gp.relevantpods, t = gp.podstarttime[p]:tstep:gp.podendtime[p]], Bin) #if 1 ==> t is departure time of pod 
	@variable(model, w_m[p = gp.relevantpods, t = gp.podstarttime[p]:tstep:gp.podendtime[p]], Bin) #if 1 ==> pod sits in queue at time t
	@variable(model, m_incomplete, Bin)

	#Objective
	@objective(model, Max, sum(sum(sum(h_m[i, p, t] for t in gp.podstarttime[p]:tstep:gp.podendtime[p]) for p in intersect(gp.relevantpods, currsol_greedysets.podswith_greedy[i])) for i in itemson[gp.m]) - 0.1 * m_incomplete - 0.001 * sum(sum(y_m[p,t] + z_m[p,t] + w_m[p,t] for t in gp.podstarttime[p]:tstep:gp.podendtime[p]) for p in gp.relevantpods) )

	#Order delivery
	if gp.orderendtime >= horizon
		@constraint(model, orderdelivery_max[i = itemson[gp.m]], sum(sum(h_m[i, p, t] for t in gp.podstarttime[p]:tstep:gp.podendtime[p]) for p in intersect(gp.relevantpods, currsol_greedysets.podswith_greedy[i])) <= 1)
		@constraint(model, orderdelivery_min[i = itemson[gp.m]], sum(sum(h_m[i, p, t] for t in gp.podstarttime[p]:tstep:gp.podendtime[p]) for p in intersect(gp.relevantpods, currsol_greedysets.podswith_greedy[i])) >= 1 - m_incomplete)
	else
		@constraint(model, orderdelivery_max[i = itemson[gp.m]], sum(sum(h_m[i, p, t] for t in gp.podstarttime[p]:tstep:gp.podendtime[p]) for p in intersect(gp.relevantpods, currsol_greedysets.podswith_greedy[i])) == 1)
	end

	#Item-pod link
	#@constraint(model, podonlyifitem[i = itemson[m], p = intersect(sp_relevantpods, podswith_greedy[i]), t = podstarttime[p]:tstep:podendtime[p], t2 = t:tstep:podendtime[p]], y_m[p,t] <= h_m[i,p,t2])
	@constraint(model, podarrives[i = itemson[gp.m], p = intersect(gp.relevantpods, currsol_greedysets.podswith_greedy[i]), t = gp.podstarttime[p]:tstep:gp.podendtime[p]], sum(y_m[p,t2] for t2 in gp.podstarttime[p]:tstep:t) >= h_m[i,p,t])
	@constraint(model, poddeparts[i = itemson[gp.m], p = intersect(gp.relevantpods, currsol_greedysets.podswith_greedy[i]), t = gp.podstarttime[p]:tstep:gp.podendtime[p]], z_m[p,t] >= h_m[i,p,t])
	@constraint(model, podwaits[p = gp.relevantpods, t = gp.podstarttime[p]:tstep:gp.podendtime[p]], w_m[p,t] >= sum(y_m[p,t2] for t2 in gp.podstarttime[p]:tstep:t) - sum(z_m[p,t2] for t2 in gp.podstarttime[p]:tstep:t))

	#Workstation throughput constraints
	@constraint(model, stationthroughput[t = sp_times_reg], itemprocesstime * sum(sum(h_m[i,p,t] for p in intersect(gp.relevantpods, currsol_greedysets.podswith_greedy[i]) if gp.podstarttime[p] <= t <= gp.podendtime[p]) for i in itemson[gp.m]) + podprocesstime * sum(z_m[p,t] for p in gp.relevantpods_t[t]) <= tstep - itemprocesstime * length(currsol.itempodpicklist[gp.w,t]) - podprocesstime * length(currsol.podsworkedat[gp.w,t]))

	#Congestion
    remainingcongestionspace = intersectiontimemaxpods - sum(currcong[p] for p in pods)
    intersectionindices = [maps.mapintersectiontorow[i] for i in currpartition.intersections] 
    timeindices = [maps.maptimetocolumn[t] for t in max(0,gp.orderstarttime):congestiontstep:min(horizon,gp.orderendtime)]
    @constraint(model, maxcongestion[i in intersectionindices, t_con in timeindices],
        sum(sum(congestionsignature[arcs[nodes[podstorageloc[p], t-arclength[gp.w,podstorageloc[p]]], nodes[gp.w,t]]][i,t_con] * y_m[p,t] for t in max(0+arclength[podstorageloc[p],gp.w],gp.podstarttime[p]):tstep:min(horizon, gp.podendtime[p])) for p in gp.relevantpods)
        + sum(sum(congestionsignature[arcs[nodes[gp.w,t], nodes[podstorageloc[p], t+arclength[gp.w,podstorageloc[p]]]]][i,t_con] * z_m[p,t] for t in max(0,gp.podstarttime[p]):tstep:min(horizon-arclength[gp.w,podstorageloc[p]], gp.podendtime[p])) for p in gp.relevantpods)
        <= max(1e-4, remainingcongestionspace[i,t_con]))

	#====================================================#

	#Solve IP
	status = optimize!(model)

	if termination_status(model) == MOI.OPTIMAL
		solvetime = solve_time(model)
		obj = 0
		for i in itemson[gp.m], p in intersect(gp.relevantpods, currsol_greedysets.podswith_greedy[i]), t in gp.podstarttime[p]:tstep:gp.podendtime[p]
			obj += value(h_m[i, p, t])
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
		for i in itemson[gp.m], p in intersect(gp.relevantpods, currsol_greedysets.podswith_greedy[i]), t in gp.podstarttime[p]:tstep:gp.podendtime[p]
			obj += value(h_m[i, p, t])
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

	return obj, solvetime, h_m, y_m, z_m, feasibleflag

end

#-----------------------------------------------------------------------------------#

function getpodstarttimes(relevantpods, orderstarttime, orderendtime, currsol, ws)

	podstarttime, podendtime = Dict(), Dict()
	sp_relevantpods = []
	sp_relevantpods_t = Dict()
	for t in orderstarttime:congestiontstep:orderendtime
		sp_relevantpods_t[t] = []
	end

	for p in relevantpods
		podtostationtime = arclength[podstorageloc[p], ws]
		podfreetimes = setdiff(orderstarttime-podtostationtime:tstep:orderendtime+podtostationtime, [t for (w,t) in currsol.podbusy[p]])
		if (podfreetimes != []) && (length(podfreetimes) == (last(podfreetimes) - podfreetimes[1]) / tstep + 1)
			push!(sp_relevantpods, p)
			if podfreetimes[1] == dummystarttime
				podstarttime[p] = max(podfreetimes[1] + podtostationtime, 0)
			else
				podstarttime[p] = max(podfreetimes[1] + podtostationtime, orderstarttime)
			end
			#### CHANGE ####
			if last(podfreetimes) >= dummyendtime
				podendtime[p] = min(last(podfreetimes) - podtostationtime, horizon)
			else
				podendtime[p] = min(last(podfreetimes) - podtostationtime, orderendtime)
			end
            for t2 in podstarttime[p]:congestiontstep:podendtime[p]
				push!(sp_relevantpods_t[t2], p)
			end
		end
	end

	return sp_relevantpods, sp_relevantpods_t, podstarttime, podendtime

end

#-----------------------------------------------------------------------------------#

#=greedyorders=[m for m in currpartition.orders if length(itemson[m]) >= 4]
greedyorders = currpartition.orders
greedypods=currpartition.pods
greedypodswith=currpartition.podswith
greedyworkstations=currpartition.workstations=#

function findgreedysolution(currpartition, currsol, greedyorders, greedypods, greedypodswith, greedyworkstations)

    #Initialize sets needed to implement greedy algorithm efficiently
    currsol_greedysets = initializeblanksolution(greedypodswith, greedypods)
    
    #Sort orders 
    sortedorders = prioritizelargerorders(greedyorders)

    #Assign orders one at a time
	totalobj = 0
	println("======= GREEDY SOLUTION =======")
	assignedorders = 0
    for m in sortedorders
		println("----- ORDER $m (", length(itemson[m]), " items) -----")
		
		#Prioritize the workstations for assignment
        overlapopportunities, workstationoverlap, relevantpods = checkpodoverlap(m, greedyworkstations, currsol_greedysets, currsol)
        totalworkstationdistance = orderclosestworkstations(m, greedyworkstations, currsol_greedysets)
        sortedworkstationtups = sort!([(w, workstationoverlap[w], totalworkstationdistance[w]) for w in greedyworkstations], by = x -> 1000*x[2] - 0.001*x[3], rev=true)
        sortedworkstations = [item[1] for item in sortedworkstationtups]
        #sortedworkstations = greedyworkstations[randperm(2)] 

		capacitybuffer2 = max(itemprocesstime, 5*min(floor(length(itemson[m])/2),3)) 
        
        #Iterate over the workstations to find an assignment
        for w2 in sortedworkstations
			
            #Find a good time block for the order-workstation combo
            feasible_flag, orderstarttime, orderendtime = fitorderin(m, w2, overlapopportunities, capacitybuffer2, currsol)
            
            if feasible_flag == 1

                #Format times and find pods with relevant items that aren't busy during the timeblock
                sp_relevantpods, sp_relevantpods_t, podstarttime, podendtime = getpodstarttimes(relevantpods, orderstarttime, orderendtime, currsol, w2)
	
                #Group together greedy problem elements into greedy problem (gp) object
                gp = (m=m, w=w2, orderstarttime=orderstarttime, orderendtime=orderendtime, relevantpods=sp_relevantpods, relevantpods_t=sp_relevantpods_t, podstarttime=podstarttime, podendtime=podendtime)
                
                #Find feasible order assignment for order m in gp
                assign_obj, assign_solvetime, h_assign, y_assign, z_assign, assignmentsuccessful = findfeasibleorderassignment(gp, currpartition, currsol, currsol_greedysets)
                
                #If good assignment, update currsol and add picks to objective
                if assignmentsuccessful == 1
					currsol, currsol_greedysets = updategreedysolution(currsol, currsol_greedysets, gp, currpartition, h_assign, y_assign, z_assign)
                    totalobj += floor(assign_obj)
					assignedorders += 1
					if visualizationflag == 1
						workstationviz(string(visualizationfolder, "/station_partition1_initial", assignedorders,"_order", m,".png"), currpartition, currsol)
					end
					break
				end
			end
		end

        if debugmode == 1
            checksolution(currpartition, currsol, debugprintstatements)
        end

    end

	return currsol

end
