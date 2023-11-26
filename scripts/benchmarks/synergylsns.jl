
#-----------------------------------------------------------------------------#

#Riley Lenaway additions

function selectingworkstations(workstationgrouplist, numsubprobworkstations, itempodpicklist, timinglist, preventcycles, tstep, horizonend, horizonstart)

	#picking first station: (Could change the order in which they are evaluated)
	#First entry: throughput in the interval (if interval is from 120 to 180, we sum throughput at 150 and 180)
	#Second Entry: orders working in the interval (same idea as above)
	#Third Entry: Total throughput combining all workstations in the time interval
	#Fourth Entry: Total throughput for the workstation over all time
	workstation_time_status = Dict()
	for i in 1:length(workstationgrouplist)
		if !(i in preventcycles)
			items_through = sum(length(itempodpicklist[workstationgrouplist[i][1],t]) for t in (timinglist[i][1]+tstep):tstep:timinglist[i][2])
			orders_open = sum(length(ordersopen[workstationgrouplist[i][1],t]) for t in (timinglist[i][1]+tstep):tstep:timinglist[i][2])
			time_throughput = sum(length(itempodpicklist[w,t]) for w in workstations for t in (timinglist[i][1]+tstep):tstep:timinglist[i][2])
			ws_throughput = sum(length(itempodpicklist[workstationgrouplist[i][1],t]) for t in 0:tstep:(horizonend-horizonstart))

			workstation_time_status[workstationgrouplist[i][1],timinglist[i][1],timinglist[i][2]] = [items_through,orders_open,time_throughput,ws_throughput]
		end
	end

	#Now we minimize the elements in the dictionary 1 by 1. So we filter the data to only include minimum of items_through,
	#then using that filtered list, we filter again to only include the min of orders_open, etc.
	x=[]
	for i in 1:length(workstationgrouplist)
		if !(i in preventcycles)
			push!(x,i)
		end
	end

	for val in 1:4
		current = 1000000
		y = deepcopy(x)
		for i in y
			value = workstation_time_status[workstationgrouplist[i][1],timinglist[i][1],timinglist[i][2]][val]
			if value == current
				push!(x,i)
			elseif value < current
				empty!(x)
				push!(x,i)
				current = value
			end
		end
	end
	#Any ties are broken randomly:
	choice = rand(x)
	sp_workstations = []
	sp_workstations = copy(workstationgrouplist[choice])
	window = timinglist[choice]

	#Now choose second workstation (do similar elimination to above for workstations during the time window):
	if numsubprobworkstations > 1
		workstation2 = Dict()
		x=[]
		for i in workstations
			if i != sp_workstations[1]
				if (i,window[1],window[2]) in keys(workstation_time_status)
					workstation2[i] = workstation_time_status[i,window[1],window[2]]
					push!(x,i)
				end
			end
		end

		for val in 1:4
			current = 1000000
			y = deepcopy(x)
			for i in y
				value = workstation2[i][val]
				if value == current
					push!(x,i)
				elseif value < current
					empty!(x)
					push!(x,i)
					current = value
				end
			end
		end

		if !(x==Any[])
			choice2 = rand(x)
			push!(sp_workstations, choice2)
		end
	end

	for i in 1:(length(preventcycles)-1)
		preventcycles[i] = preventcycles[i+1]
	end
	preventcycles[length(preventcycles)] = choice

	return sp_workstations, window#, preventcycles
end

#-----------------------------------------------------------------------------#
#sp_workstations, sp_tstart, sp_tend, tstep, remaininginventory, podstartnode, loccoords, unassignedorders, itemson, throughputperstation, targetnumpods, targetnumorders = sp_workstations, sp_tstart, sp_tend, tstep, remaininginventory, podstartnode, loccoords, unassignedorders, itemson, (podprocesstime+itemprocesstime)/60, targetnumpods, targetnumorders
function selectingordersandpods(sp_workstations, sp_tstart, sp_tend, tstep, remaininginventory, podstartnode, loccoords, assignedorders, unassignedorders, itemson, throughputperstation, targetnumpods, targetnumorders)

	#Creating a distance matrix from workstations to pods
	work_pod_dist = zeros(Float64, size(sp_workstations)[1], size(podstartnode)[1])
	for i = 1:size(work_pod_dist)[1]
		for j = 1:size(work_pod_dist)[2]
			podstart = nodelookup[podstartnode[j]][1]
			work_pod_dist[i,j] = abs(loccoords[workstations[i],1] - loccoords[podstart,1]) + abs(loccoords[workstations[i],2] - loccoords[podstart,2])
		end
	end

	#Now selecting pods
	sp_pods = []
	podscurrentlyvisiting = Dict()
	#Ensure that we include all pods already visiting the workstations at times within the neighborhood:
	for w in sp_workstations
		order_list = []
		for t in max(0,sp_tstart):tstep:max(horizon,sp_tend)
			for i in 1:length(itempodpicklist[w,t])
				if !(itempodpicklist[w,t][i][3] in sp_pods)
					push!(sp_pods,itempodpicklist[w,t][i][3])
					push!(order_list,itempodpicklist[w,t][i][3])
				end
			end
		end
		podscurrentlyvisiting[w] = deepcopy(order_list)
	end

	#Now select some pods closest (within 10 meters) to each workstation (pod-workstation synergy)
	toadd = (targetnumpods - length(sp_pods))/2
	counter = 0
	if size(sp_workstations)[1] > 1
		for i = 1:size(work_pod_dist)[2]
			if counter%6 == 0 && round(counter/6) <= toadd
				if work_pod_dist[1,i] < 10 || work_pod_dist[2,i] < 10
					if !(i in sp_pods)
						push!(sp_pods,i)
					end
				end
			end
			counter += 1
		end
	end
	if size(sp_workstations)[1] == 1
		for i = 1:size(work_pod_dist)[2]
			if counter%3 == 0 && round(counter/3) <= toadd
				if work_pod_dist[1,i] < 10
					if !(i in sp_pods)
						push!(sp_pods,i)
					end
				end
			end
			counter += 1
		end
	end

	#Now we select orders to add to the neighborhood
	sp_orders = assignedorders
	numordersremaining = targetnumorders - length(assignedorders)

	#First, for each workstation, select orders with the highest % of items already on pods visiting (order-workstation synergy)
	#(max 1 per "size" per station)
	smallorderpodpercent = Dict()
	mediumorderpodpercent = Dict()
	largeorderpodpercent = Dict()
	for w in sp_workstations
		for j in unassignedorders
			order_perc = []
			for k in podscurrentlyvisiting[w]
				count=0
				for h in itemson[j]
					if (h,k) in keys(remaininginventory)
						count=count+1
					end
				end
				push!(order_perc, count/length(itemson[j]))
			end
			if order_perc != []
				if length(itemson[j]) <= 2
					smallorderpodpercent[w,j] = maximum(order_perc)
				elseif length(itemson[j]) >2 && length(itemson[j]) < 8
					mediumorderpodpercent[w,j] = maximum(order_perc)
				else
					largeorderpodpercent[w,j] = maximum(order_perc)
				end
			end
		end
	end

	#(We now add 1*# of workstations orders for each order size)
	for i in 1:length(sp_workstations)
		if !(podscurrentlyvisiting == Dict()) && length(sp_orders) != length(unassignedorders)
			if length(smallorderpodpercent) != 0
				maxsmall = findmax(smallorderpodpercent)
				sp_orders = union(sp_orders, maxsmall[2][2])
				smallorderpodpercent[maxsmall[2]] = 0
			end
			if length(mediumorderpodpercent) != 0
				maxmedium = findmax(mediumorderpodpercent)
				sp_orders = union(sp_orders, maxmedium[2][2])
				mediumorderpodpercent[maxmedium[2]] = 0
			end
			if length(largeorderpodpercent) != 0
				maxlarge = findmax(largeorderpodpercent)
				sp_orders = union(sp_orders, maxlarge[2][2])
				largeorderpodpercent[maxlarge[2]] = 0
			end
		end
	end

	#Next, select orders that have the highest percentage of their items on the average pod in our subproblem (order-pod synergy)
	smallorderpodpercent = Dict()
	mediumorderpodpercent = Dict()
	largeorderpodpercent = Dict()
	for j in unassignedorders
		if !(j in sp_orders)
			order_perc = []
			for k in pods
				count=0
				for h in itemson[j]
					if (h,k) in keys(remaininginventory)
						count=count+1
					end
				end
				push!(order_perc, count/length(itemson[j]))
			end
			if order_perc != []
				if length(itemson[j]) <= 2
					smallorderpodpercent[j] = mean(order_perc)
				elseif length(itemson[j]) >2 && length(itemson[j]) <6
					mediumorderpodpercent[j] = mean(order_perc)
				else
					largeorderpodpercent[j] = mean(order_perc)
				end
			end
		end
	end

	#We now add orders until we reach our goal for the subproblem
	i=0
	while i == 0
		if length(sp_orders) < length(unassignedorders) && length(sp_orders) < targetnumorders
			if (length(smallorderpodpercent) != 0) & (length(sp_orders) < targetnumorders)
				maxsmall = findmax(smallorderpodpercent)
				push!(sp_orders, maxsmall[2])
				smallorderpodpercent[maxsmall[2]] = 0
			end
			if (length(mediumorderpodpercent) != 0) & (length(sp_orders) < targetnumorders)
				maxmedium = findmax(mediumorderpodpercent)
				push!(sp_orders, maxmedium[2])
				mediumorderpodpercent[maxmedium[2]] = 0
			end
			if (length(largeorderpodpercent) != 0) & (length(sp_orders) < targetnumorders)
				maxlarge = findmax(largeorderpodpercent)
				push!(sp_orders, maxlarge[2])
				largeorderpodpercent[maxlarge[2]] = 0
			end
		else
			i=1
		end
	end

	#Now add more pods with most items from the subproblem orders (more order pod synergy)
	all_sp_items = []
	for i in sp_orders, j in itemson[i]
		push!(all_sp_items, j)
	end
	podscore = Dict()
	for p in pods
		if !(p in sp_pods)
			count = 0
			for i in all_sp_items
				if (i,p) in keys(remaininginventory)
					count=count+1
				end
			end
			podscore[p] = count
		end
	end

	#Limit to max number of pods
	while length(sp_pods) < targetnumpods
		push!(sp_pods, findmax(podscore)[2])
		podscore[findmax(podscore)[2]] = 0
	end

	return sp_pods, sp_orders

end

#-----------------------------------------------------------------------------#

function writefinalsolutionoutputs(outputfilename, obj_full, solvetime_init, solvetime_sp, solvetime_build, solvetime_update, h_currsol, y_currsol, v_currsol)

	itemspicked = 0
	itemspickedfrom = Dict()
	podpicklist = Dict()
	for m in orders
		itemspickedfrom[m] = []
	end
	for w in workstations, t in times
		podpicklist[w,t] = []
	end
	for w in workstations, t in times
		for (m,i,p) in itempodpicklist[w,t]
			itemspicked += 1
			push!(itemspickedfrom[m], i)
			podpicklist[w,t] = union(podpicklist[w,t], p)
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
		l = min(11,convert(Int64, length(itemson[m])))
		orderfreq[l] += 1
	end

	order_open_time_per_item = sum(sum(sum(v_currsol[m,w,t] for t in times) for m in orderscompleted) for w in workstations) / sum(length(itemson[m]) for m in orderscompleted)

	congestionat = Dict()
	for l in intersections, t in 0:congestiontstep:horizon
		congestionat[l,t] = 0
		for p in pods, a in intersect(congestionarcs[l,t], podarcset[p])
			congestionat[l,t] += congestioncontribution[a,l,t] * y_currsol[p,a]
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
			method = [string(methodname, "_", orderprioritization)], 
			objective = [obj_full], 
			solve_time = [solvetime_init + solvetime_sp + solvetime_build + solvetime_update],
			solvetime_init = [solvetime_init],
			solvetime_sp = [solvetime_sp],
			solvetime_build = [solvetime_build],
			solvetime_update = [solvetime_update],
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

#-----------------------------------------------------------------------------#

function writegreedyoutputs(outputfilename, obj_full, solvetime_greedy, h_currsol, y_currsol, v_currsol)

	itemspicked = 0
	itemspickedfrom = Dict()
	podpicklist = Dict()
	for m in orders
		itemspickedfrom[m] = []
	end
	for w in workstations, t in times
		podpicklist[w,t] = []
	end
	for w in workstations, t in times
		for (m,i,p) in itempodpicklist[w,t]
			itemspicked += 1
			push!(itemspickedfrom[m], i)
			podpicklist[w,t] = union(podpicklist[w,t], p)
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
		l = min(11,convert(Int64, length(itemson[m])))
		orderfreq[l] += 1
	end

	order_open_time_per_item = sum(sum(sum(v_currsol[m,w,t] for t in times) for m in orderscompleted) for w in workstations) / sum(length(itemson[m]) for m in orderscompleted)

	congestionat = Dict()
	for l in intersections, t in 0:congestiontstep:horizon
		congestionat[l,t] = 0
		for p in pods, a in intersect(congestionarcs[l,t], podarcset[p])
			congestionat[l,t] += congestioncontribution[a,l,t] * y_currsol[p,a]
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
			method = [string(methodname, "_", orderprioritization)], 
			objective = [obj_full], 
			solve_time = [solvetime_greedy ],
			solvetime_greedy = [solvetime_greedy],
			solvetime_null1 = [0],
			solvetime_null2 = [0],
			solvetime_null3 = [0],
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


#-----------------------------------------------------------------------------#

function writelsnsoutput(outputfilename, objlist)

	df = DataFrame(run_id = [run_id for j in 1:length(objlist)],
			test_instance_id = [instance_id for j in 1:length(objlist)],
			warehouse_id = [warehouse_id for j in 1:length(objlist)], 
			random_seed = [random_seed for j in 1:length(objlist)], 
			method = [string(methodname, "_", orderprioritization) for j in 1:length(objlist)], 
			lsnsiteration = [j for j in 1:length(objlist)],
			objvalue = objlist
           )

	CSV.write(outputfilename, df)

end