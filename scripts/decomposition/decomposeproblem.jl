
function createproblempartitions(stationsperpartition)

	stationsin, storagelocsin, intersectionsin, podsin = Dict(), Dict(), Dict(), Dict()

	numpartitions = convert(Int, ceil(numworkstations/stationsperpartition))
	partitions = 1:numpartitions

	for s in partitions
		stationsin[s] = []
		storagelocsin[s] = []
		intersectionsin[s] = []
		podsin[s] = []
	end

	if warehouse_y_length_meters > 45
		numpartitions_x = numpartitions/2
		partition_x_length = warehouse_x_length_meters / numpartitions_x
		partition_y_length = warehouse_y_length_meters / 2
	else
		numpartitions_x = numpartitions
		partition_x_length = warehouse_x_length_meters / numpartitions_x
		partition_y_length = warehouse_y_length_meters 
	end	

	#Station and location distribution
	partitionindex = 1
	for x in 1:numpartitions_x, y in 1:2
		min_x, max_x = (x-1) * partition_x_length, x * partition_x_length - 1e-4
		min_y, max_y = (y-1) * partition_y_length, y * partition_y_length - 1e-4
		#println("Partition $partitionindex is x = ($min_x, $max_x) and y = ($min_y, $max_y)")

		for w in workstations
			if (min_x <= loccoords[w,1] <= max_x) & (min_y <= loccoords[w,2] <= max_y)
				push!(stationsin[partitionindex], w)
			end
		end

		for l in storagelocs
			if (min_x <= loccoords[l,1] <= max_x) & (min_y <= loccoords[l,2] <= max_y)
				push!(storagelocsin[partitionindex], l)
			end
		end

		for l in intersections
			if (min_x <= intcoords[l][1] <= max_x) & (min_y <= intcoords[l][2] <= max_y)
				push!(intersectionsin[partitionindex], l)
			end
		end

		partitionindex +=1
	end

	#Pod distribution
	for s in partitions, l in storagelocsin[s], p in intersect(pods,podsstoredat[l]) 
		push!(podsin[s], p)
	end

	return numpartitions, partitions, stationsin, storagelocsin, intersectionsin, podsin

end

#----------------------------------------------------------------------------------------------------#

function maximizesynergyoforderassignments(partitionobjective, beta, features, featureinfo, featurenums, numpartitions, partitions, stationsin, podsin)

	totalitems = sum(length(itemson[m]) for m in orders)
	itembuffer = totalitems * 0.08
	allpartitions = 1:numpartitions+1

	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", 60) # iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", 0)
	set_optimizer_attribute(model, "MIPGap", 0.01)

	#Variables
	@variable(model, x[orders, allpartitions], Bin)
	@variable(model, itembuff >= 0)

	#Objective = maximize synergy 
	if partitionobjective == "random"
		@objective(model, Max, sum(sum(rand()*x[m,s] for m in orders) for s in allpartitions))
	elseif partitionobjective == "synergy"
		println("Doing God's work")
		@objective(model, Max, sum(sum(sum(sum(beta.mp[f] * features.mp[chop(featureinfo.names[f],head=3,tail=0)][m,p] * x[m,s] for f in featurenums.mp) for m in orders) for p in podsin[s]) for s in partitions)) 
	elseif partitionobjective == "heuristic"
		@objective(model, Min, 2*itembuff + sum(x[m,numpartitions+1] for m in orders)) #sum(sum(sum(sum(beta_mp[f] * mp_features[chop(featurenames[f],head=3,tail=0)][m,p] * x[m,s] for f in mp_featnums) for m in orders) for p in podsin[s]) for s in partitions) )
	elseif partitionobjective == "none"
		@objective(model, Max, 0)
		@constraint(model, [m in orders], x[m,1] == 1)
	end

	#Constraints
	@constraint(model, orderassignment[m in orders], sum(x[m,s] for s in allpartitions) == 1)
	@constraint(model, possibilityoffulfillment[m in orders, i in itemson[m], s in partitions], x[m,s] <= length(intersect(podswith[i], podsin[s])))
	@constraint(model, partitionbalance[s in partitions], sum(length(itemson[m]) * x[m,s] for m in orders) >= totalitems/numpartitions - itembuff) 
	
	#====================================================#

	#Solve IP
	status = optimize!(model)

	#====================================================#

	if (termination_status(model) == MOI.OPTIMAL) || ((termination_status(model) == MOI.TIME_LIMIT) && (has_values(model)))
		return value.(x), solve_time(model)
	elseif (termination_status(model) == MOI.INFEASIBLE) #|| (termination_status(model) == MOI.INF_OR_UNBD)
		println("Infeasible SP selection!")
		return None, 0
	else
		println("No solution found!")
		println(termination_status(model))
		return None, 0
	end

end

#----------------------------------------------------------------------------------------------------#

function distributeordersbysynergy(x, partitions, numpartitions)

	ordersin = Dict()
	itemsin = []

	for s in partitions
		ordersin[s] = []
		for m in orders
			if value(x[m,s]) > 1e-4
				push!(ordersin[s], m)
			end
		end
	end

	dummypartition = numpartitions + 1 
	unassignedorders = []
	for m in orders
		if value(x[m,dummypartition]) > 1e-4
			push!(unassignedorders, m)
		end
	end

	for s in partitions
		push!(itemsin, [])
		for m in ordersin[s]
			itemsin[s] = union(itemsin[s], itemson[m])
		end
	end

	return ordersin, unassignedorders, itemsin

end

#----------------------------------------------------------------------------------------------------#

function additionalsets(podsin, partitions)

	podswith_s = []
	allitemsin = []
	for s in partitions
		newdict = Dict()
		for p in podsin[s], i in podstartinventory[p]
			try
				newdict[i] = union(newdict[i], p)
			catch
				newdict[i] = [p]
			end
		end
		push!(podswith_s, newdict)
		push!(allitemsin, keys(newdict))
	end

	return podswith_s, allitemsin
	
end


#----------------------------------------------------------------------------------------------------#

function decomposeproblem(stationsperpartition, partitionobjective, beta, features, featureinfo, featurenums)

	numpartitions, partitions, stationsin, storagelocsin, intersectionsin, podsin = createproblempartitions(stationsperpartition)
	orderassignments, partitionsolvetime = maximizesynergyoforderassignments(partitionobjective, beta, features, featureinfo, featurenums, numpartitions, partitions, stationsin, podsin)
	ordersin, unassignedorders, itemsin = distributeordersbysynergy(orderassignments, partitions, numpartitions)
	podswith_s, allitemsin = additionalsets(podsin, partitions)

	partitioninfo = []
	for s in partitions
		part = (partitionid=s, workstations=stationsin[s], storagelocs=storagelocsin[s], intersections=intersectionsin[s], pods=podsin[s], orders=ordersin[s], podswith=podswith_s[s], allitems=allitemsin[s], items=itemsin[s])
		push!(partitioninfo, part)
	end
	
	globalpartitionid = numpartitions+1
	part = (partitionid=globalpartitionid, workstations=workstations, storagelocs=storagelocs, 
	intersections=intersectionsin, pods=pods, orders=unassignedorders, podswith=podswith, 
	allitems=allitems, items=items)
	push!(partitioninfo, part)

	return numpartitions, partitions, partitioninfo, globalpartitionid, partitionsolvetime

end

