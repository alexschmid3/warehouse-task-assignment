
function writeglobalsolutionoutputs(globalsolutionfilename, solvemetrics)

	objective, full_solve_time, solvetime_init, solvetime_sp, time_utilization, throughput_utilization, bestthroughput_utilization, total_orders_worked, total_orders_completed, congestion_utilization, max_congestion, pods_used, unique_pods_used, items_picked_per_pod, pod_distance_travelled, order_open_time_per_item, orders_size_1, orders_size_2, orders_size_3, orders_size_4, orders_size_5, orders_size_6, orders_size_7, orders_size_8, orders_size_9, orders_size_10, orders_size_large = zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1)
	totalpodtrips, totalpodstops, multistoptrips, multistoppods, multitrippods, usefulitemsperpod = zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1), zeros(numpartitions+1)

	usefulnumerator, usefuldenominator = 0, 0
	for s in 1:numpartitions
		currpartition = partitioninfo[s]
		currsol = partitionsolution[s]
		objective[s] += sum(sum(length(currsol.itempodpicklist[w,t]) for t in times) for w in currpartition.workstations)
		time_utilization[s] += (podprocesstime * sum(sum(length(currsol.podsworkedat[w,t]) for t in times) for w in currpartition.workstations) + itemprocesstime * sum(sum(length(currsol.itempodpicklist[w,t]) for t in times) for w in currpartition.workstations)) / ((tstep+horizon) * length(currpartition.workstations))
		throughput_utilization[s] += ((podprocesstime + itemprocesstime) * sum(sum(length(currsol.itempodpicklist[w,t]) for t in times) for w in currpartition.workstations)) / (horizon * length(currpartition.workstations))
		bestthroughput_utilization[s] += (podprocesstime * sum(sum(min(1,length(currsol.podsworkedat[w,t])) for t in times) for w in currpartition.workstations) + itemprocesstime * sum(sum(length(currsol.itempodpicklist[w,t]) for t in times) for w in currpartition.workstations)) / (horizon * length(currpartition.workstations))

		ordersworked, podsworked, itemsdone = [], [], Dict()
		for w in currpartition.workstations, t in times, (m,i,p) in currsol.itempodpicklist[w,t]
			ordersworked = union(ordersworked, m)
			try
				itemsdone[m] += 1
			catch
				itemsdone[m] = 1
			end
		end
		
		for w in currpartition.workstations, t in times, p in currsol.podsworkedat[w,t]
			push!(podsworked, p)
			#totalpoddist += 2*warehousedistance[podstorageloc[p],w]
		end
	
		totalpoddist = 0
		for p in currpartition.pods
			for a in currsol.ypath[p]
				l1,l2 = nodelookup[arclookup[a][1]][1], nodelookup[arclookup[a][2]][1]
				totalpoddist += warehousedistance[l1,l2]
			end
		end
		orderscompleted = [m for m in ordersworked if itemsdone[m] >= length(itemson[m])]
		total_orders_worked[s] += length(ordersworked)
		total_orders_completed[s] += length(orderscompleted)

		intindices, timeindices = [maps.mapintersectiontorow[i] for i in currpartition.intersections], [maps.maptimetocolumn[t] for t in 0:congestiontstep:horizon]
		allcongestionvalues = sum(currcong[p][intindices, timeindices] for p in currpartition.pods)
		allcongestionmaxes = map(Int, (intersectionmaxpods[l] for l=currpartition.intersections, j=1:length(0:congestiontstep:horizon))) 
		congestionutilization = allcongestionvalues ./ allcongestionmaxes
		congestion_utilization[s] += sum(allcongestionvalues) / sum(allcongestionmaxes)
		max_congestion[s] += maximum(congestionutilization) 

		pods_used[s] += length(podsworked)
		unique_pods_used[s] += length(unique!(podsworked))
		items_picked_per_pod[s] = objective[s] / pods_used[s]
		pod_distance_travelled[s] += totalpoddist

		orderopentime = Dict()
		for m in currpartition.orders
			orderopentime[m] = 1
		end
		for w in currpartition.workstations, t in times, m in currsol.ordersopen[w,t] 
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
	
		#Reporting for multi-stop analysis
		totalstopsperpod_dict = Dict()
		totaltripsperpod_dict = Dict()
		for p in currpartition.pods
			totalstopsperpod_dict[p] = 0
			totaltripsperpod_dict[p] = sum(sum(partitionsolution[s].y[p,a] for a in intersect(A_plus[extendednodes[podstorageloc[p],t]], A_space, podarcset[p])) for t in -tstep:tstep:horizon) 
		end
		for t in times, w in currpartition.workstations, p in unique([p[3] for p in [item for item in partitionsolution[s].itempodpicklist[w,t]]]) #p in partitionsolution[s].podsworkedat[w,t]
			totalstopsperpod_dict[p] += 1
		end
		for p in currpartition.pods
			totaltripsperpod_dict[p] = min(totaltripsperpod_dict[p], totalstopsperpod_dict[p])
		end
		totalstopsperpod = collect(values(totalstopsperpod_dict))
		totaltripsperpod = collect(values(totaltripsperpod_dict)) #[sum(sum(partitionsolution[s].y[p,a] for a in intersect(A_plus[nodes[podstorageloc[p],t]], A_space)) for t in times) for p in currpartition.pods]

		totalpodtrips[s] += sum(totaltripsperpod)
		totalpodstops[s] += sum(totalstopsperpod) 
		multistoptrips[s] += totalpodstops[s] - totalpodtrips[s]
		multistoppods[s] += length([p for p in currpartition.pods if totalstopsperpod_dict[p] >= 2])
		multitrippods[s] += length([p for p in currpartition.pods if totaltripsperpod_dict[p] >= 2])

		#Useful items
		usefulitems = [length(intersect(items, podstartinventory[p])) for p in currpartition.pods]
		usefulitemsperpod[s] += sum(usefulitems) / sum(num_items_per_pod for p in currpartition.pods)
		usefulnumerator += sum(usefulitems) 
		usefuldenominator += sum(num_items_per_pod for p in currpartition.pods)

	end
	usefulitemsperpod[numpartitions+1] += usefulnumerator / usefuldenominator

	df = DataFrame(
			row_id = [row_id for s in 1:numpartitions+1],
			run_id = [run_id for s in 1:numpartitions+1],
			test_instance_id = [instance_id for s in 1:numpartitions+1],
			warehouse_id = [warehouse_id for s in 1:numpartitions+1],
			random_seed = [random_seed for s in 1:numpartitions+1],
			method = [methodname for s in 1:numpartitions+1],
			partition = union([string(s) for s in 1:numpartitions], ["global"]),
			lsnsiteration = push!(["summary" for p in 1:numpartitions], string(20)),
			objective = objective,
			full_solve_time = solvemetrics.solve_time,
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
			unique_pods_used = unique_pods_used,
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
			orders_size_large = orders_size_large,
			totalpodtrips = totalpodtrips, 
			totalpodstops = totalpodstops, 
			multistoptrips = multistoptrips, 
			multistoppods = multistoppods, 
			multitrippods = multitrippods,
			usefulitemsperpod = usefulitemsperpod
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
			full_solve_time = [],
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
			orders_size_large = [],
			totalpodtrips = [], 
			totalpodstops = [], 
			multistoptrips = [], 
			multistoppods = [], 
			multitrippods = [],
			usefulitemsperpod = []
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
	allcongestionmaxes = map(Int, (intersectionmaxpods[l] for l=currpartition.intersections, j=1:length(0:congestiontstep:horizon))) 
	congestionutilization = allcongestionvalues ./ allcongestionmaxes	
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
			full_solve_time = [iterationtime],
			solvetime_init = [inittime],
			solvetime_spselection = [spselectiontime],
			solvetime_spopt = [sp_solvetime],
			time_utilization = [(podprocesstime * sum(sum(length(currsol.podsworkedat[w,t]) for t in times) for w in currpartition.workstations) + itemprocesstime * sum(sum(length(currsol.itempodpicklist[w,t]) for t in times) for w in currpartition.workstations)) / ((tstep+horizon) * length(currpartition.workstations))],
			throughput_utilization = [((podprocesstime + itemprocesstime) * sum(sum(length(currsol.itempodpicklist[w,t]) for t in times) for w in currpartition.workstations)) / (horizon * length(currpartition.workstations))],
			bestthroughput_utilization = [(podprocesstime * sum(sum(min(1,length(currsol.podsworkedat[w,t])) for t in times) for w in currpartition.workstations) + itemprocesstime * sum(sum(length(currsol.itempodpicklist[w,t]) for t in times) for w in currpartition.workstations)) / (horizon * length(currpartition.workstations))],
			total_orders_worked = [length(ordersworked)],
			total_orders_completed = [length(orderscompleted)],
			congestion_utilization = [ sum(allcongestionvalues) / sum(allcongestionmaxes)],
			max_congestion = [maximum(congestionutilization)],
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
			orders_size_large = length([m for m in orderscompleted if length(itemson[m]) >= 11]),
			totalpodtrips = [0], 
			totalpodstops = [0], 
			multistoptrips = [0], 
			multistoppods = [0], 
			multitrippods = [0],
			usefulitemsperpod = [0]
		)	

	CSV.write(globalsolutionfilename, df, append=true)

end


#-----------------------------------------------------------------------------------#

function writeglobalsolutionoutputs_partitioning(partitionsolvetime)

	df = DataFrame(
			row_id = [row_id],
			run_id = [run_id],
			test_instance_id = [instance_id],
			warehouse_id = [warehouse_id],
			random_seed = [random_seed],
			method = [methodname],
			partition = ["global"],
			lsnsiteration = ["partitioningproblem"],
			objective = [0],
			full_solve_time = [partitionsolvetime],
			solvetime_init = [0],
			solvetime_spselection = [0],
			solvetime_spopt = [0],
			time_utilization = [0],
			throughput_utilization = [0],
			bestthroughput_utilization = [0],
			total_orders_worked = [0],
			total_orders_completed = [0],
			congestion_utilization = [0],
			max_congestion = [0],
			pods_used = [0],
			unique_pods_used  = [0],
			items_picked_per_pod = [0],
			pod_distance_travelled = [0],
			order_open_time_per_item = [0],
			orders_size_1 = [0],
			orders_size_2 = [0],
			orders_size_3 = [0],
			orders_size_4 = [0],
			orders_size_5 = [0],
			orders_size_6 = [0],
			orders_size_7 = [0],
			orders_size_8 = [0],
			orders_size_9 = [0],
			orders_size_10 = [0],
			orders_size_large = [0],
			totalpodtrips = [0], 
			totalpodstops = [0], 
			multistoptrips = [0], 
			multistoppods = [0], 
			multitrippods = [0],
			usefulitemsperpod = [0]
		)	

	CSV.write(globalsolutionfilename, df, append=true)

end

