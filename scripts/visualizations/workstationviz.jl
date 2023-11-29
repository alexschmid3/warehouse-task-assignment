
using Luxor, Colors

tempHueList = collect(Colors.color_names)

#Exclude colors that are too light
hueList = []
for item in tempHueList
	if (item[2][1] + item[2][2] + item[2][3] <= 515) & (max(abs(item[2][1]-item[2][2]), abs(item[2][1]-item[2][3]), abs(item[2][2]-item[2][3])) >= 100)
		push!(hueList, item) 
	end
end

gradientHueList = []
r1,g1,b1 = 139,209,249 
r2,g2,b2 = 0,114,178
#r1,g1,b1 = 237,125,49 
#r2,g2,b2 = 46,117,182
push!(gradientHueList, (237,125,49))
for j in 1:10
	rj = convert(Int64, round(r1 - (j - 1) * (r1-r2) / 9, digits=0))
	gj = convert(Int64, round(g1 - (j - 1) * (g1-g2) / 9, digits=0))
	bj = convert(Int64, round(b1 - (j - 1) * (b1-b2) / 9, digits=0))
	push!(gradientHueList, (rj, gj, bj))
end

colorindex = 1
ordersizecolor = Dict()
for m in orders
	if length(itemson[m]) == 1
		ordersizecolor[m] = gradientHueList[1]
	else
		ordersizecolor[m] = gradientHueList[min(11, 1+length(itemson[m]))]
	end
	global colorindex += 1
end

#---------------------------------------------------------------------------------------#

function workstationviz(wsdrawingname, partition, currsol)

	#Find coordinates for each time-space node
	w_adj = minimum(partition.workstations) - 2
	nodePoints = Dict()
	slotPoints = Dict()
	numtimesteps = horizon/tstep + 1
	timesperrow = ceil(numtimesteps / 2)
	k1 = 1000/(timesperrow + 1)  
	k2 = 600/(length(partition.workstations) + 3)

	for w in partition.workstations, t in 0:tstep:horizon
		n = nodes[w,t]
		ycoord = (w - w_adj) + (floor((t/tstep)/timesperrow))/2
		xcoord = mod((t/tstep), timesperrow) + 1 + floor((t/tstep)/timesperrow)
		tup = (-400 + xcoord*k1 - 65, -350 + ycoord*k2)   #-250, -350
		nodePoints[n] = Point(tup)
		slotPoints[n] = Dict()
	end

	#---------------------------------------------------------------------------------------#

	throughput = Dict()
	for w in partition.workstations, t in 0:tstep:horizon
		sortedlist = sort(currsol.itempodpicklist[w,t], by = x -> x[3])
		currentpod = -1
		vizlist = []
		for item in sortedlist
			if item[3] == currentpod
				push!(vizlist, item[1])
			else
				push!(vizlist, "P")
				push!(vizlist, item[1])
				currentpod = item[3]
			end
		end
		throughput[w,t] = vizlist
	end

	#---------------------------------------------------------------------------------------#

	#Create arc list including arc properties (start point, end point, color, dash) 
	#arcList, labelList = [], []
	#outputCounter = Dict()
	#queueCounter = Dict()
	#for w in workstations, t in times
	#	queueCounter[nodes[w,t]] = 0
	#end

	#Order segments
	#for sgmt in podsegments	
	#	p = sgmt[1]
	#	if (sgmt[4] in workstations) & (sgmt[5] in storagelocs)
	#		startPoint = nodePoints[nodes[sgmt[4], sgmt[2]]]
	#		endPoint = nodePoints[nodes[sgmt[4], sgmt[2]+tstep]]
	#		try
	#			outputCounter[(startPoint, endPoint)] += 1
	#		catch
	#			outputCounter[(startPoint, endPoint)] = 1
	#		end
	#	elseif (sgmt[4] in storagelocs) & (sgmt[5] in workstations)
	#		nd = nodePoints[nodes[sgmt[5], sgmt[2]]]
	#		try
	#			queueCounter[nd] += 1
	#		catch
	#			queueCounter[nd] = 1
	#		end
	#	end
	#end

	#for item in outputCounter
	#	startPoint, endPoint = item[1]
	#	newColor = "black"
	#	newDash = "solid"
	#	push!(arcList, (startPoint, endPoint, newColor, outputCounter[(startPoint, endPoint)], newDash))
	#end

	#Create throughput labels 
	#for w in workstations, t in tstep:tstep:horizon-tstep
	#	if solutionthroughput[w,t] > 0.01
	#		n1, n2 = nodes[w,t], nodes[w,t+tstep]
	#		nodePoints[n1]
	#		startPoint, endPoint = nodePoints[n1], nodePoints[n2]
	#		newColor = "black"
	#		newDash = "solid"
	#		push!(arcList, (startPoint, endPoint, newColor, convert(Int64,solutionthroughput[w,t]), newDash))
	#		push!(labelList, (startPoint, endPoint, convert(Int64,solutionthroughput[w,t])))
	#	end
	#end

	#-------------------------------------------------------------------------#

	#orderslots = []
	#or w in workstations, t in 0:tstep:horizon
	#	n = nodes[w,t]
	#	local y = nodelookup[n][1] - numstoragelocs + 1
	#	local x = (nodelookup[n][2]/tstep)+1
	#	for os in 1:workstationordercapacity
	#		tup = (-400 + x*k1, -350 + y*k2 + os*min(30,k2/(workstationordercapacity+2)) + 5) 
	#		startPoint = Point(tup[1] - 6, tup[2]) 
	#		endPoint = Point(tup[1] + 6, tup[2]) 
	#		linecolor = "black"
	#		linedash = "solid"
	#		push!(orderslots, (startPoint, endPoint, linecolor, linedash))
	#	end
	#end

	#-------------------------------------------------------------------------#

	#Order-station assignment
	orderlistbynode = Dict()
	orderslotassign = Dict()
	slotorderassign = Dict()
	for w in partition.workstations, t in times, s in 1:workstationordercapacity
		slotorderassign[nodes[w,t], s] = 0 
	end
	for w in partition.workstations, t in times
		orderlistbynode[w,t] = currsol.ordersopen[w,t]
	end
	for w in partition.workstations, t in times
		n = nodes[w,t]
		orderlist = orderlistbynode[w,t]
		if t >= tstep
			for m in intersect(orderlist, orderlistbynode[w,t-tstep])
				try
					s = orderslotassign[nodes[w,t-tstep],m] 
					orderslotassign[n,m] = s
					slotorderassign[n,s] = m
				catch
					for s in 1:workstationordercapacity
						if slotorderassign[n,s] == 0
							orderslotassign[n,m] = s
							slotorderassign[n,s] = m
							break
						end
					end
				end
			end
			for m in setdiff(orderlist, orderlistbynode[w,t-tstep])
				for s in 1:workstationordercapacity
					if slotorderassign[n,s] == 0
						orderslotassign[n,m] = s
						slotorderassign[n,s] = m
						break
					end
				end
			end
		else
			for m in orderlist
				for s in 1:workstationordercapacity
					if slotorderassign[n,s] == 0
						orderslotassign[n,m] = s
						slotorderassign[n,s] = m
						break
					end
				end
			end
		end
	end

	#Orders at stations
	orderdots = []
	for item in orderslotassign
		n, m, s = item[1][1], item[1][2], item[2]
		local y = nodelookup[n][1] - numstoragelocs + 1
		local x = (nodelookup[n][2]/tstep)+1 
		dotpoint = Point((-400 + x*k1, -350 + y*k2 + s*min(30,k2/(workstationordercapacity+2))))
		dotcolor = hueList[1]#mod(m+14,533)][2]
		push!(orderdots, (dotpoint, dotcolor))
	end

	#-------------------------------------------------------------------------#

	Drawing(1200,700, wsdrawingname)
	origin()
	background("white")

	#Draw arcs
	#fontsize(12)
	#for i in arcList
	#	setcolor(i[3])
	#	setdash(i[5])
	#	line(i[1], i[2] , :stroke)
		
	#	theta = atan((i[2][2] - i[1][2])/(i[2][1] - i[1][1]))
	#	dist = distance(i[1], i[2])
		
	#	arrowhead = (1-8/dist)*i[2] + (8/dist)*i[1] #8 pixels from the end node
		
	#	local p = ngon(arrowhead, 5, 3, theta, vertices=true)
	#	poly(p, :fill,  close=true)
	#end

	#Nodes
	setcolor("black")
	setdash("solid")
	lineseparation = (nodePoints[nodes[partition.workstations[1],tstep]] - nodePoints[nodes[partition.workstations[1],0]])[1] / 6
	linelength = 0.8 * lineseparation
	setline(2)
	for w in partition.workstations, t in times
		np = nodePoints[nodes[w,t]]
		#Luxor.circle(np, 4, :fill)
		for t2 in 0:5
			np1 = Point(np[1] + lineseparation * t2, np[2] + 10)
			np2 = Point(np[1] + lineseparation * t2 + linelength, np[2] + 10)
			Luxor.line(np1, np2, :stroke)
			slotPoints[nodes[w,t]][t2+1] = Point(np[1] + lineseparation * t2 + linelength / 2, np[2] + 5)
		end
	end

	#Add subproblem highlights
	#setline(7)
	#setcolor("goldenrod1")
	#for w in sp_workstations, t in intersect(0:tstep:horizon, sp_times)
	#	np = nodePoints[nodes[w,t]]
	#	np1 = Point(np[1] , np[2] + 16)
	#	np2 = Point(np[1] + k1, np[2] + 16)   
	#	Luxor.line(np1, np2, :stroke)
	#end

	fontsize(16)
	setcolor("black")
	#Location labels
	for l in partition.workstations
		coord = nodePoints[nodes[(l,0.0)]]
		Luxor.label("Station $l     ", :W , coord)
	end

	fontsize(12)
	#Time labels
	setline(1)
	for w in partition.workstations, t in times
		coord = nodePoints[nodes[(w,t)]] + Point(0,-k2/10)
		if isinteger(t/60)
			setcolor("black")
			Luxor.label("t = $t", :N , coord)
		end
		setcolor("gray75")
		Luxor.line(nodePoints[nodes[(w,t)]] + Point(0,11), nodePoints[nodes[(w,t)]] - Point(0,11), :stroke)
	end

	fontsize(8)
	for w in partition.workstations, t in times
		for j in 1:length(throughput[w,t])
			if throughput[w,t][j] == "P"
				setline(0.5)
				setcolor("gray35")
				#Luxor.text("P", slotPoints[nodes[w,t]][j])
				Luxor.squircle(slotPoints[nodes[w,t]][j], 3, 3, rt=0.0, :stroke)
				Luxor.line(slotPoints[nodes[w,t]][j] + Point(-3,-3), slotPoints[nodes[w,t]][j] + Point(3,3), :stroke)
				Luxor.line(slotPoints[nodes[w,t]][j] + Point(-3,3), slotPoints[nodes[w,t]][j] + Point(3,-3), :stroke)
			else
				#setcolor("black")
				r_val, g_val, b_val = ordersizecolor[throughput[w,t][j][1]]
				setcolor(convert(Colors.HSV, Colors.RGB(r_val / 255, g_val / 255, b_val / 255)))
				#mycolor = ordercolor[throughput[w,t][j][1]]
				Luxor.squircle(slotPoints[nodes[w,t]][j], 3, 3, rt=0.0, :fillpreserve)
			end
		end
	end

	#Arc labels
	#setcolor("black")
	#for lbl in labelList
	#	coord = 0.5 * (lbl[1] + lbl[2]) + Point(0,-3)
	#	label(string(lbl[3]), :N , coord)
	#end

	#Order slots
	#for os in orderslots
	#	setcolor(os[3])
	#	setdash(os[4])
	#	line(os[1], os[2] , :stroke)
	#end

	#Orders in slots
	#for d in orderdots
	#	r_val, g_val, b_val = d[2][1]/255, d[2][2]/255, d[2][3]/255
	#	setcolor(convert(Colors.HSV, Colors.RGB(r_val, g_val, b_val)))
	#	circle(d[1], 3, :fill)
	#end

	#Queue labels
	#setcolor("gray50")
	#for w in workstations, t in times
	#	coord = nodePoints[nodes[(w,t)]] + Point(0,max(-k2/6, -15))
	#	if queueCounter[nodes[(w,t)]] > 0.01
	#		label(string("Q = ", queueCounter[nodes[(w,t)]]), :N , coord)
	#	end
	#end

	finish()
	preview()
	
	return nothing

end

#---------------------------------------------------------------------------------------#

function workstationviz_three(wsdrawingname, itempodpicklist)

	#Find coordinates for each time-space node
	nodePoints = Dict()
	slotPoints = Dict()
	numtimesteps = horizon/tstep + 1
	timesperrow = ceil(numtimesteps / 3)
	k1 = 1000/(timesperrow + 1)  
	k2 = 600/(numworkstations + 2)

	for w in workstations, t in 0:tstep:horizon
		n = nodes[w,t]
		if nodelookup[n][2]/tstep < timesperrow
			ymod = (floor((nodelookup[n][2]/tstep)/timesperrow))/3
		else
			ymod = (floor((nodelookup[n][2]/tstep - 1)/(timesperrow-1)))/3
		end 
		ycoord = (nodelookup[n][1] - numstoragelocs + 1) + ymod
		if nodelookup[n][2]/tstep < timesperrow
			xmod = mod((nodelookup[n][2]/tstep), timesperrow)
		else
			xmod = mod((nodelookup[n][2]/tstep)-1, timesperrow-1)+1
		end 
		xcoord = xmod + 1 #+ floor((nodelookup[n][2]/tstep)/timesperrow)
		tup = (-400 + xcoord*k1 - 65, -450 + ycoord*k2)   #-250, -350
		nodePoints[n] = Point(tup)
		slotPoints[n] = Dict()
	end

	#-------------------------------------------------------------------------#

	throughput = Dict()
	for w in workstations, t in 0:tstep:horizon
		sortedlist = sort(itempodpicklist[w,t], by = x -> x[3])
		currentpod = -1
		vizlist = []
		for item in sortedlist
			if item[3] == currentpod
				push!(vizlist, item[1])
			else
				push!(vizlist, "P")
				push!(vizlist, item[1])
				currentpod = item[3]
			end
		end
		throughput[w,t] = vizlist
	end

	#-------------------------------------------------------------------------#

	#Order-station assignment
	orderlistbynode = Dict()
	orderslotassign = Dict()
	slotorderassign = Dict()
	for w in workstations, t in times, s in 1:workstationordercapacity
		slotorderassign[nodes[w,t], s] = 0 
	end
	for n in 1:numnodes
		try
			orderlistbynode[n] = openorders[n]
		catch
			orderlistbynode[n] = []
		end
	end
	for n in 1:numnodes
		orderlist = orderlistbynode[n]
		w, t = nodelookup[n]
		if t >= tstep
			for m in intersect(orderlist, orderlistbynode[nodes[w,t-tstep]])
				try
					s = orderslotassign[nodes[w,t-tstep],m] 
					orderslotassign[n,m] = s
					slotorderassign[n,s] = m
				catch
					for s in 1:workstationordercapacity
						if slotorderassign[n,s] == 0
							orderslotassign[n,m] = s
							slotorderassign[n,s] = m
							break
						end
					end
				end
			end
			for m in setdiff(orderlist, orderlistbynode[nodes[w,t-tstep]])
				for s in 1:workstationordercapacity
					if slotorderassign[n,s] == 0
						orderslotassign[n,m] = s
						slotorderassign[n,s] = m
						break
					end
				end
			end
		else
			for m in orderlist
				for s in 1:workstationordercapacity
					if slotorderassign[n,s] == 0
						orderslotassign[n,m] = s
						slotorderassign[n,s] = m
						break
					end
				end
			end
		end
	end

	#Orders at stations
	orderdots = []
	for item in orderslotassign
		n, m, s = item[1][1], item[1][2], item[2]
		local y = nodelookup[n][1] - numstoragelocs + 1
		local x = (nodelookup[n][2]/tstep)+1 
		dotpoint = Point((-400 + x*k1, -350 + y*k2 + s*min(30,k2/(workstationordercapacity+2))))
		dotcolor = hueList[mod(m+14,533)][2]
		push!(orderdots, (dotpoint, dotcolor))
	end

	#-------------------------------------------------------------------------#

	Drawing(1200,700, wsdrawingname)
	origin()
	background("white")

	#Nodes
	setcolor("black")
	setdash("solid")
	lineseparation = (nodePoints[nodes[workstations[1],tstep]] - nodePoints[nodes[workstations[1],0]])[1] / 6
	linelength = 0.8 * lineseparation
	setline(2)
	for w in workstations, t in times
		np = nodePoints[nodes[w,t]]
		#Luxor.circle(np, 4, :fill)
		for t2 in 0:5
			np1 = Point(np[1] + lineseparation * t2, np[2] + 10)
			np2 = Point(np[1] + lineseparation * t2 + linelength, np[2] + 10)
			Luxor.line(np1, np2, :stroke)
			slotPoints[nodes[w,t]][t2+1] = Point(np[1] + lineseparation * t2 + linelength / 2, np[2] + 3)
		end
	end

	fontsize(16)
	setcolor("black")
	#Location labels
	for l in workstations
		coord = nodePoints[nodes[(l,0.0)]]
		Luxor.label("Station $l     ", :W , coord)
	end

	fontsize(12)
	#Time labels
	setline(1)
	for w in workstations, t in times
		coord = nodePoints[nodes[(w,t)]] + Point(0,-k2/10)
		if t == 0		
			setcolor("black")
			Luxor.label("t = $t", :N , coord)
		elseif isinteger(t/60)
			setcolor("black")
			Luxor.label("$t", :N , coord)
		end
		setcolor("gray75")
		Luxor.line(nodePoints[nodes[(w,t)]] + Point(0,11), nodePoints[nodes[(w,t)]] - Point(0,11), :stroke)
	end

	fontsize(8)
	for w in workstations, t in times
		for j in 1:length(throughput[w,t])
			if throughput[w,t][j] == "P"
				setline(0.5)
				setcolor("gray55")
				#Luxor.text("P", slotPoints[nodes[w,t]][j])
				Luxor.squircle(slotPoints[nodes[w,t]][j], 4, 4, rt=0.0, :stroke)
				Luxor.line(slotPoints[nodes[w,t]][j] + Point(-4,-4), slotPoints[nodes[w,t]][j] + Point(4,4), :stroke)
				Luxor.line(slotPoints[nodes[w,t]][j] + Point(-4,4), slotPoints[nodes[w,t]][j] + Point(4,-4), :stroke)
			else
				#setcolor("black")
				r_val, g_val, b_val = ordersizecolor[throughput[w,t][j][1]]
				setcolor(convert(Colors.HSV, Colors.RGB(r_val / 255, g_val / 255, b_val / 255)))
				#mycolor = ordercolor[throughput[w,t][j][1]]
				#setcolor(mycolor[1])
				Luxor.squircle(slotPoints[nodes[w,t]][j], 5, 5, rt=0.0, :fillpreserve)
			end
		end
	end

	finish()
	preview()
	
	return nothing

end

