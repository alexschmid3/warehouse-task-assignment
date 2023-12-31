
function writeglobalsolutionoutputs(globalsolutionfilename, solvemetrics)

	objective, solve_time, solvetime_init, solvetime_sp, time_utilization, throughput_utilization, bestthroughput_utilization, total_orders_worked, total_orders_completed, congestion_utilization, max_congestion, pods_used, unique_pods_used, items_picked_per_pod, pod_distance_travelled, order_open_time_per_item, orders_size_1, orders_size_2, orders_size_3, orders_size_4, orders_size_5, orders_size_6, orders_size_7, orders_size_8, orders_size_9, orders_size_10, orders_size_large = zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1)

	for s in 1:numpartitions
		currpartition = partitioninfo[s]
		objective[s] += sum(sum(length(partitionsolution[s].itempodpicklist[w,t]) for t in times) for w in currpartition.workstations)
		time_utilization[s] += (podprocesstime * sum(sum(length(partitionsolution[s].podsworkedat[w,t]) for t in times) for w in currpartition.workstations) + itemprocesstime * sum(sum(length(partitionsolution[s].itempodpicklist[w,t]) for t in times) for w in currpartition.workstations)) / ((tstep+horizon) * length(currpartition.workstations))
		throughput_utilization[s] += ((podprocesstime + itemprocesstime) * sum(sum(length(partitionsolution[s].itempodpicklist[w,t]) for t in times) for w in currpartition.workstations)) / (horizon * length(currpartition.workstations))
		bestthroughput_utilization[s] += (podprocesstime * sum(sum(min(1,length(partitionsolution[s].podsworkedat[w,t])) for t in times) for w in currpartition.workstations) + itemprocesstime * sum(sum(length(partitionsolution[s].itempodpicklist[w,t]) for t in times) for w in currpartition.workstations)) / (horizon * length(currpartition.workstations))

		ordersworked, podsworked, itemsdone = [], [], Dict()
		for w in currpartition.workstations, t in times, (m,i,p) in partitionsolution[s].itempodpicklist[w,t]
			ordersworked = union(ordersworked, m)
			try
				itemsdone[m] += 1
			catch
				itemsdone[m] = 1
			end
		end
		totalpoddist = 0
		for w in currpartition.workstations, t in times, p in partitionsolution[s].podsworkedat[w,t]
			push!(podsworked, p)
			totalpoddist += 2*warehousedistance[podstorageloc[p],w]
		end
		orderscompleted = [m for m in ordersworked if itemsdone[m] >= length(itemson[m])]
		total_orders_worked[s] += length(ordersworked)
		total_orders_completed[s] += length(orderscompleted)

		intindices, timeindices = [maps.mapintersectiontorow[i] for i in currpartition.intersections], [maps.maptimetocolumn[t] for t in 0:congestiontstep:horizon]
		allcongestionvalues = sum(currcong[p][intindices, timeindices] for p in currpartition.pods)
		congestion_utilization[s] += mean(allcongestionvalues)
		max_congestion[s] += maximum(allcongestionvalues) 

		pods_used[s] += length(podsworked)
		unique_pods_used[s] += length(unique!(podsworked))
		items_picked_per_pod[s] = objective[s] / pods_used[s]
		pod_distance_travelled[s] += totalpoddist

		orderopentime = Dict()
		for m in ordersworked
			orderopentime[m] = 1
		end
		for w in currpartition.workstations, t in times, m in partitionsolution[s].ordersopen[w,t] 
			orderopentime[m] += 1
		end
		if length(ordersworked) > 0
			order_open_time_per_item[s] += sum([orderopentime[m] for m in ordersworked]) / sum([length(itemson[m]) for m in ordersworked])
		else
			order_open_time_per_item[s] = 0
		end

		orders_size_1[s] = length([m for m in orderscompleted if length(itemson[m]) == 1])
		orders_size_2[s] = length([m for m in orderscompleted if length(itemson[m]) == 2])
		orders_size_3[s] = length([m for m in orderscompleted if length(itemson[m]) == 3])
		orders_size_4[s] = length([m for m in orderscompleted if length(itemson[m]) == 4])
		orders_size_5[s] = length([m for m in orderscompleted if length(itemson[m]) == 5])
		orders_size_6[s] = length([m for m in orderscompleted if length(itemson[m]) == 6])
		orders_size_7[s] = length([m for m in orderscompleted if length(itemson[m]) == 7])
		orders_size_8[s] = length([m for m in orderscompleted if length(itemson[m]) == 8])
		orders_size_9[s] = length([m for m in orderscompleted if length(itemson[m]) == 9])
		orders_size_10[s] = length([m for m in orderscompleted if length(itemson[m]) == 10])
		orders_size_large[s] = length([m for m in orderscompleted if length(itemson[m]) >= 11])

	end

	df = DataFrame(
			row_id = [row_id for s in 1:numpartitions+1],
			run_id = [run_id for s in 1:numpartitions+1],
			test_instance_id = [instance_id for s in 1:numpartitions+1],
			warehouse_id = [warehouse_id for s in 1:numpartitions+1],
			random_seed = [random_seed for s in 1:numpartitions+1],
			method = [methodname for s in 1:numpartitions+1],
			partition = union([string(s) for s in 1:numpartitions], ["global"]),
			lsnsiteration = union(["summary" for p in 1:numpartitions], [string(20)]),
			objective = objective,
			solve_time = solvemetrics.solve_time,
			solvetime_init = solvemetrics.solvetime_init,
			solvetime_spselection = solvemetrics.solvetime_spsel,
			solvetime_spopt = solvemetrics.solvetime_sp,
			time_utilization = time_utilization,
			throughput_utilization = throughput_utilization,
			bestthroughput_utilization = bestthroughput_utilization,
			total_orders_worked = total_orders_worked,
			total_orders_completed = total_orders_completed,
			congestion_utilization = congestion_utilization,
			max_congestion = max_congestion,
			pods_used = pods_used,
			unique_pods_used  = unique_pods_used,
			items_picked_per_pod = items_picked_per_pod,
			pod_distance_travelled = pod_distance_travelled,
			order_open_time_per_item = order_open_time_per_item,
			orders_size_1 = orders_size_1,
			orders_size_2 = orders_size_2,
			orders_size_3 = orders_size_3,
			orders_size_4 = orders_size_4,
			orders_size_5 = orders_size_5,
			orders_size_6 = orders_size_6,
			orders_size_7 = orders_size_7,
			orders_size_8 = orders_size_8,
			orders_size_9 = orders_size_9,
			orders_size_10 = orders_size_10,
			orders_size_large = orders_size_large
		)

	CSV.write(globalsolutionfilename, df, append=true)

end

#-----------------------------------------------------------------------------------#

function writeglobalsolutionoutputs_init(globalsolutionfilename)

	df = DataFrame(
			row_id = [],
			run_id = [],
			test_instance_id = [],
			warehouse_id = [],
			random_seed = [],
			method = [],
			partition = [],
			lsnsiteration = [],
			objective = [],
			solve_time = [],
			solvetime_init = [],
			solvetime_spselection = [],
			solvetime_spopt = [],
			time_utilization = [],
			throughput_utilization = [],
			bestthroughput_utilization = [],
			total_orders_worked = [],
			total_orders_completed = [],
			congestion_utilization = [],
			max_congestion = [],
			pods_used = [],
			unique_pods_used  = [],
			items_picked_per_pod = [],
			pod_distance_travelled = [],
			order_open_time_per_item = [],
			orders_size_1 = [],
			orders_size_2 = [],
			orders_size_3 = [],
			orders_size_4 = [],
			orders_size_5 = [],
			orders_size_6 = [],
			orders_size_7 = [],
			orders_size_8 = [],
			orders_size_9 = [],
			orders_size_10 = [],
			orders_size_large = []
		)

	CSV.write(globalsolutionfilename, df)

end

#-----------------------------------------------------------------------------------#

function writeglobalsolutionoutputs_iter(sp_iter, inittime, iterationtime, spselectiontime, sp_solvetime, globalsolutionfilename, s, currpartition, currsol)

	ordersworked, podsworked, itemsdone = [], [], Dict()
	for w in currpartition.workstations, t in times, (m,i,p) in currsol.itempodpicklist[w,t]
		ordersworked = union(ordersworked, m)
		try
			itemsdone[m] += 1
		catch
			itemsdone[m] = 1
		end
	end
	totalpoddist = 0
	for w in currpartition.workstations, t in times, p in currsol.podsworkedat[w,t]
		push!(podsworked, p)
		totalpoddist += 2*warehousedistance[podstorageloc[p],w]
	end
	orderscompleted = [m for m in ordersworked if itemsdone[m] >= length(itemson[m])]
	intindices, timeindices = [maps.mapintersectiontorow[i] for i in currpartition.intersections], [maps.maptimetocolumn[t] for t in 0:congestiontstep:horizon]
	allcongestionvalues = sum(currcong[p][intindices, timeindices] for p in currpartition.pods)
	orderopentime = Dict()
	for m in currpartition.orders
		orderopentime[m] = 0
	end
	for m in ordersworked
		orderopentime[m] += 1
	end
	for w in currpartition.workstations, t in times, m in currsol.ordersopen[w,t] 
		orderopentime[m] += 1
	end
	if length(ordersworked) > 0
		order_open_time_per_item = sum([orderopentime[m] for m in ordersworked]) / sum([length(itemson[m]) for m in ordersworked])
	else
		order_open_time_per_item = 0
	end
	objective = sum(sum(length(currsol.itempodpicklist[w,t]) for t in times) for w in currpartition.workstations)

	df = DataFrame(
			row_id = [row_id],
			run_id = [run_id],
			test_instance_id = [instance_id],
			warehouse_id = [warehouse_id],
			random_seed = [random_seed],
			method = [methodname],
			partition = [string(s)],
			lsnsiteration = [sp_iter],
			objective = [objective],
			solve_time = [iterationtime],
			solvetime_init = [inittime],
			solvetime_spselection = [spselectiontime],
			solvetime_spopt = [sp_solvetime],
			time_utilization = [(podprocesstime * sum(sum(length(currsol.podsworkedat[w,t]) for t in times) for w in currpartition.workstations) + itemprocesstime * sum(sum(length(currsol.itempodpicklist[w,t]) for t in times) for w in currpartition.workstations)) / ((tstep+horizon) * length(currpartition.workstations))],
			throughput_utilization = [((podprocesstime + itemprocesstime) * sum(sum(length(currsol.itempodpicklist[w,t]) for t in times) for w in currpartition.workstations)) / (horizon * length(currpartition.workstations))],
			bestthroughput_utilization = [(podprocesstime * sum(sum(min(1,length(currsol.podsworkedat[w,t])) for t in times) for w in currpartition.workstations) + itemprocesstime * sum(sum(length(currsol.itempodpicklist[w,t]) for t in times) for w in currpartition.workstations)) / (horizon * length(currpartition.workstations))],
			total_orders_worked = [length(ordersworked)],
			total_orders_completed = [length(orderscompleted)],
			congestion_utilization = [mean(allcongestionvalues)],
			max_congestion = [maximum(allcongestionvalues)],
			pods_used = [length(podsworked)],
			unique_pods_used  = [length(unique!(podsworked))],
			items_picked_per_pod = [objective / length(podsworked)],
			pod_distance_travelled = [totalpoddist],
			order_open_time_per_item = [order_open_time_per_item],
			orders_size_1 = length([m for m in orderscompleted if length(itemson[m]) == 1]),
			orders_size_2 = length([m for m in orderscompleted if length(itemson[m]) == 2]),
			orders_size_3 = length([m for m in orderscompleted if length(itemson[m]) == 3]),
			orders_size_4 = length([m for m in orderscompleted if length(itemson[m]) == 4]),
			orders_size_5 = length([m for m in orderscompleted if length(itemson[m]) == 5]),
			orders_size_6 = length([m for m in orderscompleted if length(itemson[m]) == 6]),
			orders_size_7 = length([m for m in orderscompleted if length(itemson[m]) == 7]),
			orders_size_8 = length([m for m in orderscompleted if length(itemson[m]) == 8]),
			orders_size_9 = length([m for m in orderscompleted if length(itemson[m]) == 9]),
			orders_size_10 = length([m for m in orderscompleted if length(itemson[m]) == 10]),
			orders_size_large = length([m for m in orderscompleted if length(itemson[m]) >= 11])
		)	

	CSV.write(globalsolutionfilename, df, append=true)

end

