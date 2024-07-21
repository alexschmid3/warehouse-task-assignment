
include("findassignedorders.jl")

#-----------------------------------------------------------------------------------------------------#

function selectingworkstations(currpartition, windows, currsol, numsubprobworkstations, tabulist)

	#picking first station: (Could change the order in which they are evaluated)
	#First entry: throughput in the interval (if interval is from 120 to 180, we sum throughput at 150 and 180)
	#Second Entry: orders working in the interval (same idea as above)
	#Third Entry: Total throughput combining all workstations in the time interval
	#Fourth Entry: Total throughput for the workstation over all time
	workstation_time_status = Dict()
	x=[]
	for i in 1:length(windows)
		win = windows[i]
		if !(i in tabulist) & (length(win.workstations) == numsubprobworkstations)
			items_through = sum(sum(length(currsol.itempodpicklist[w,t]) for w in win.workstations) for t in win.times)
			orders_open = sum(sum(length(currsol.ordersopen[w,t]) for w in win.workstations) for t in win.times)
			time_throughput = sum(sum(length(currsol.itempodpicklist[w,t]) for w in currpartition.workstations) for t in win.times)
			ws_throughput = sum(sum(length(currsol.itempodpicklist[w,t]) for w in win.workstations) for t in 0:tstep:horizon)

			workstation_time_status[i] = (items_through, orders_open, time_throughput, ws_throughput)
			push!(x,i)
		end
	end

	#Now we minimize the elements in the dictionary 1 by 1. So we filter the data to only include minimum of items_through,
	#then using that filtered list, we filter again to only include the min of orders_open, etc.
	
	for val in 1:4
		current = 1000000
		y = deepcopy(x)
		for i in y
			value = workstation_time_status[i][val]
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
	sp_window = windows[choice]

	#Update tabu list
	if length(tabulist) >= maxtabu
		popfirst!(tabulist)
	end
	push!(tabulist, choice)

	return choice, sp_window, tabulist
end

#-----------------------------------------------------------------------------------------------------#

function selectingordersandpods(sp_window, currpartition, currsol, targetnumpods, targetnumorders)

	#Creating a distance matrix from workstations to pods
	work_pod_dist = Dict() # zeros(Float64, size(sp_window.workstations)[1], size(podstartnode)[1])
	for i in sp_window.workstations, j in currpartition.pods
		podstart = nodelookup[podstartnode[j]][1]
		work_pod_dist[i,j] = abs(loccoords[i,1] - loccoords[podstart,1]) + abs(loccoords[i,2] - loccoords[podstart,2])
	end

	#Now selecting pods
	sp_pods = []
	podscurrentlyvisiting = Dict()
	#Ensure that we include all pods already visiting the workstations at times within the neighborhood:
	for w in sp_window.workstations
		order_list = []
		for t in sp_window.times
			for i in 1:length(currsol.itempodpicklist[w,t])
				if !(currsol.itempodpicklist[w,t][i][3] in sp_pods)
					push!(sp_pods,currsol.itempodpicklist[w,t][i][3])
					push!(order_list,currsol.itempodpicklist[w,t][i][3])
				end
			end
		end
		podscurrentlyvisiting[w] = deepcopy(order_list)
	end

	#Now select some pods closest (within 10 meters) to each workstation (pod-workstation synergy)
	toadd = (targetnumpods - length(sp_pods))/2
	counter = 0
	if size(sp_window.workstations)[1] > 1
		for i in currpartition.pods
			if counter%6 == 0 && round(counter/6) <= toadd
				if minimum([work_pod_dist[w,i] for w in sp_window.workstations]) < 10 
					if !(i in sp_pods)
						push!(sp_pods,i)
					end
				end
			end
			counter += 1
		end
	end
	if size(sp_window.workstations)[1] == 1
		for i in currpartition.pods
			if counter%3 == 0 && round(counter/3) <= toadd
				if work_pod_dist[sp_window.workstations[1],i] < 10
					if !(i in sp_pods)
						push!(sp_pods,i)
					end
				end
			end
			counter += 1
		end
	end

	#Now we select orders to add to the neighborhood
	assignedorders, unassignedorders, numadditionalorders = findassignedorders(currpartition, sp_window, currsol, targetnumorders)
	sp_orders = copy(assignedorders)

	#First, for each workstation, select orders with the highest % of items already on pods visiting (order-workstation synergy)
	#(max 1 per "size" per station)
	smallorderpodpercent = Dict()
	mediumorderpodpercent = Dict()
	largeorderpodpercent = Dict()
	for w in sp_window.workstations
		for j in unassignedorders
			order_perc = []
			for k in podscurrentlyvisiting[w]
				count=0
				for h in itemson[j]
					if (h,k) in keys(inventory)
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
	for i in 1:length(sp_window.workstations)
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
					if (h,k) in keys(inventory)
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
			if (length([m for m in values(smallorderpodpercent) if m > -0.5]) != 0) & (length(sp_orders) < targetnumorders)
				maxsmall = findmax(smallorderpodpercent)
				push!(sp_orders, maxsmall[2])
				smallorderpodpercent[maxsmall[2]] = -1
			end
			if (length([m for m in values(mediumorderpodpercent) if m > -0.5]) != 0) & (length(sp_orders) < targetnumorders)
				maxmedium = findmax(mediumorderpodpercent)
				push!(sp_orders, maxmedium[2])
				mediumorderpodpercent[maxmedium[2]] = -1
			end
			if (length([m for m in values(largeorderpodpercent) if m > -0.5]) != 0) & (length(sp_orders) < targetnumorders)
				maxlarge = findmax(largeorderpodpercent)
				push!(sp_orders, maxlarge[2])
				largeorderpodpercent[maxlarge[2]] = -1
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
	for p in currpartition.pods
		if !(p in sp_pods)
			count = 0
			for i in all_sp_items
				if (i,p) in keys(inventory)
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

	#Find the relevant items - those that are in the assigned orders AND are picked during the time window
	sp_itemson = Dict()
	for m in assignedorders
		sp_itemson[m] = []
		for i in itemson[m]
			#if (sum(sum(currsol.h[m,i,p,w,sp_window.tend] for p in currpartition.podswith[i]) for w in sp_window.workstations) > 0.01) && (sp_window.tstart >= 0) && (sum(sum(currsol.h[m,i,p,w,sp_window.tstart-tstep] for p in currpartition.podswith[i]) for w in sp_window.workstations) < 0.01)
			#	push!(sp_itemson[m], i)
			if checkiteminpicklist(m, i, max(0,sp_window.tstart), min(horizon,sp_window.tend), sp_window.workstations, currsol)
				push!(sp_itemson[m], i)
			#	println("FOUND AN ISSUE = $m, $i")
			elseif (sp_window.tstart == -1*tstep) && (sum(sum(currsol.h[m,i,p,w,sp_window.tend] for p in currpartition.podswith[i]) for w in sp_window.workstations) > 0.01)
				push!(sp_itemson[m], i)
			end
		end
	end

	additionalorders = setdiff(sp_orders, assignedorders)
	for m in additionalorders
		sp_itemson[m] = copy(itemson[m])
	end

	sp_items = []
	for m in sp_orders
		sp_items = union(sp_items, sp_itemson[m])
	end

	return sp_pods, sp_orders, sp_itemson, sp_items

end

#-----------------------------------------------------------------------------------------------------#

function selectsynergisticsubproblem(currpartition, windows, currsol, targetnumorders, targetnumpods, numsubprobworkstations, tabulist)

	sp_winid, sp_window, tabulist = selectingworkstations(currpartition, windows, currsol, numsubprobworkstations, tabulist)
	sp_pods, sp_orders, sp_itemson, sp_items = selectingordersandpods(sp_window, currpartition, currsol, targetnumpods, targetnumorders)

	return sp_winid, sp_orders, sp_window, sp_pods, sp_itemson, sp_items, tabulist

end