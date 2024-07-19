
function updatelastoptimizeddifference(lastoptimizeddifference, tabulist, win, windows, predicted_obj, actual_obj)

	for w in setdiff(1:length(windows), tabulist)
		lastoptimizeddifference[w] = 0
	end
	if lastoptpenaltyflag == 1
		lastoptimizeddifference[win] = abs(predicted_obj - actual_obj)
	end

	return lastoptimizeddifference

end

#-----------------------------------------------------------------------------------------------------#

function reoptimizesubproblem(sp, currsol, currpartition, turnoffcongestion_flag, ipoutputflag=0, iptimelimit=300)

	model = Model(() -> Gurobi.Optimizer(GRB_ENV)) #Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", timeforreooptimization) # iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", 0)  #ipoutputflag)
		
	#Variables
	buildtimestart=time()
	if (debugmode == 1) & (debugprintstatements == 1)
		println("Orders = ", sp.orders)
		println("Pods = ", sp.pods)
		println("Stations = ", sp.workstations)
		println("Times = ", sp.times)
		for m in sp.ordersinprogress
			println("-------------------------------------------------------------")
			println("Order $m - ", sp.itemson[m]," vs. " , itemson[m])
			println("  Assigned --> ", currsol.stationassign[m])
			println("-------------------------------------------------------------")
		end
		for m in sp.ordersinprogress, i in itemson[m], p in currpartition.podswith[i], w in currpartition.workstations
			for t in times
				if currsol.h[m,i,p,w,t] > 1e-4
					println("h[$m, $i, $p, $w, $t] = ", currsol.h[m,i,p,w,t])
					break
				end
			end
		end
	end

	@variable(model, h[m = sp.orders, i = sp.itemson[m], p = sp.podswith[i], sp.workstations, sp.times], Bin)
	@variable(model, y[p = sp.pods, a = sp.podarcset[p]], Bin)
	@variable(model, z[sp.orders, sp.workstations], Bin)
	@variable(model, v[sp.orders, sp.workstations, sp.times], Bin)
	@variable(model, f[sp.orders, sp.workstations, sp.times], Bin)
	@variable(model, g[sp.orders, sp.workstations, sp.times], Bin)

	#Objective
	@objective(model, Max, sum(sum(sum(sum(h[m, i, p, w, last(sp.times)] for w in sp.workstations) for p in sp.podswith[i]) for i in sp.itemson[m]) for m in sp.orders) )

	#Order assignment and delivery constraints
	asgn=time()
	@constraint(model, stationassignment[m = sp.orders], sum(z[m, w] for w in sp.workstations) <= 1)
	@constraint(model, orderdelivery_inprogress[m = sp.ordersinprogress, i = sp.itemson[m]], sum(h[m, i, p, currsol.stationassign[m], last(sp.times)] for p in sp.podswith[i]) == 1)
	@constraint(model, orderdelivery_max[m = sp.orders, i = sp.itemson[m], w = sp.workstations], sum(h[m, i, p, w, last(sp.times)] for p in sp.podswith[i]) == z[m, w])
	if debugmode == 1
		println("   Assignment constraint time = ", time()-asgn)
	end

	#Linking constraints
	lnk=time()
	@constraint(model, podlink[m = sp.orders, i = sp.itemson[m], p = sp.podswith[i], w = sp.workstations, t = setdiff(sp.times,sp.times[1])], h[m,i,p,w,t] - h[m,i,p,w,t-tstep] <= sum(y[p,a] for a in setdiff(intersect(A_plus_p[p,nodes[w,t]], sp.arcset), A_queues)) + sp.y_known[p,nodes[w,t]] )
	@constraint(model, podlinkzero[m = sp.orders, i = sp.itemson[m], p = sp.podswith[i], w = sp.workstations], h[m,i,p,w,sp.times[1]] <= sum(y[p,a] for a in setdiff(intersect(A_plus_p[p,nodes[w,sp.times[1]]], sp.arcset), A_queues)) + sp.y_known[p,nodes[w,sp.times[1]]] )
	@constraint(model, workstationlink[m = sp.orders, i = sp.itemson[m], p = sp.podswith[i], w = sp.workstations, t = setdiff(sp.times,sp.times[1])], h[m,i,p,w,t] - h[m,i,p,w,t-tstep] <= z[m, w] )
	@constraint(model, nondecreasing[m = sp.orders, i = sp.itemson[m], p = sp.podswith[i], w = sp.workstations, t = setdiff(sp.times,sp.times[1])], h[m,i,p,w,t] >= h[m,i,p,w,t-tstep] )
	if debugmode == 1
		println("   Linking constraint time = ", time()-lnk)
	end

	#Inventory constraints
	inv=time()
	@constraint(model, maxinventory[i = sp.items, p = sp.podswith[i]], sum(sum(h[m,i,p,w,last(sp.times)] for w in sp.workstations) for m in intersect(sp.orders, sp.orderswith[i])) <= sp.remaininginventory[i,p])
	if debugmode == 1
		println("   Inventory constraint time = ", time()-inv)
	end

	#Pod movement constraints
	pm=time()
	@constraint(model, podflowbalance[p = sp.pods, n in intersect(podnodeset[p], sp.nodeset)], sum(y[p,a] for a in intersect(A_minus_p[p,n], sp.arcset)) - sum(y[p,a] for a in intersect(A_plus_p[p,n], sp.arcset)) == sp.podsupply[p,n])
	if debugmode == 1
		println("   Pod flow constraint time = ", time()-pm)
	end

	#Workstation capacity constraints
	oop=time()
	@constraint(model, fconfig_init[m = sp.orders, w = sp.workstations, t = sp.times], f[m,w,t] >= currsol.v[m,w,sp.times[1]-tstep])
	@constraint(model, fconfig[m = sp.orders, i = sp.itemson[m], p = sp.podswith[i], w = sp.workstations, t = sp.times], f[m,w,t] >= h[m,i,p,w,t])
	@constraint(model, gconfig_init[m = sp.orders, i = sp.itemson[m], w = sp.workstations, t = sp.times], g[m,w,t] <= 1 - currsol.v[m,w,last(sp.times)])
	@constraint(model, gconfig[m = sp.orders, i = sp.itemson[m], w = sp.workstations, t = sp.times], g[m,w,t] <= sum(h[m,i,p,w,t] for p in sp.podswith[i]))
	@constraint(model, vconfig[m = sp.orders, w = sp.workstations, t = sp.times], v[m,w,t] == f[m,w,t] - g[m,w,t])
	@constraint(model, stationcapacity[w = sp.workstations, t = sp.times], sum(v[m,w,t] for m in sp.orders) <= C[w] )
	if debugmode == 1
		println("   Orders open constraint time = ", time()-oop)
	end
	#Workstation throughput constraints
	tpt=time()
	@constraint(model, stationthroughput[w = sp.workstations, t = setdiff(sp.times, sp.times[1])], itemprocesstime * sum(sum(sum(h[m,i,p,w,t] - h[m,i,p,w,t-tstep] for p in sp.podswith[i]) for i in sp.itemson[m]) for m in sp.orders) + podprocesstime * sum(sum(y[p,a] for a in intersect(A_plus_p[p,nodes[w,t]], A_space, sp.arcset)) for p in sp.pods) + podprocesstime * sum(sp.y_known[p, nodes[w,t]] for p in sp.pods) <= tstep)
	@constraint(model, stationthroughput_first[w = sp.workstations, t = sp.times[1]], itemprocesstime * sum(sum(sum(h[m,i,p,w,t] for p in sp.podswith[i]) for i in sp.itemson[m]) for m in sp.orders) + podprocesstime * sum(sum(y[p,a] for a in intersect(A_plus_p[p,nodes[w,t]], A_space, sp.arcset)) for p in sp.pods) + podprocesstime * sum(sp.y_known[p, nodes[w,t]] for p in sp.pods) <= tstep)
	if debugmode == 1
		println("   Throughput constraint time = ", time()-tpt)
	end

	#Congestion
	#@constraint(model, maxcongestion[l in partitionfloor.reducedintersections, t in max(0,sp.tstart):congestiontstep:min(horizon,sp.tend)], sum(sum(partitionfloor.congestioncontribution[a,l,t] * y[p,a] for a in intersect(partitionfloor.congestionarcs[l,t], sp.podarcset[p], sp.arcset)) for p in sp.pods) <= intersectionmaxpods[l] - sp.ambientcongestion[l,t])
	#@constraint(model, maxcongestion[l in currpartition.intersections, t in max(0,sp.tstart):congestiontstep:min(horizon,sp.tend)], sum(congestionsignature[a][maps.mapintersectiontorowcolumn[l], maps.maptimetocol[t]] * y[p,a] for p in sp.pods) <= intersectionmaxpods[l] - sp.ambientcongestion[l,t])
	ct1=time()
	if turnoffcongestion_flag == 0
		@constraint(model, maxcongestion[l in currpartition.intersections, t in max(0,sp.tstart):congestiontstep:min(horizon,sp.tend)], sum(sum(congestionsignature[a][maps.mapintersectiontorow[l], maps.maptimetocolumn[t]] * y[p,a] for a in sp.sp_podarcset_cong[p]) for p in sp.pods) <= intersectionmaxpods[l] - sum(currcong[p][maps.mapintersectiontorow[l],maps.maptimetocolumn[t]] for p in setdiff(pods, sp.pods)))
	end
	if debugmode == 1
		println("   Congestion constraint time = ", time()-ct1)
	end

	#Storage location capacity
	if anystoragelocation_flag == 1
		println("Adding storage loc capacity constraint...")
		@time @constraint(model, storageloccapacity[s in storagelocs, t in 0:tstep:horizon-tstep], sum(y[p,arcs[nodes[s,t],nodes[s,t+tstep]]] for p in sp.pods) + sum(currsol.y[p,arcs[nodes[s,t],nodes[s,t+tstep]]] for p in setdiff(currpartition.pods,sp.pods)) <= numpods / numstoragelocs)
	end

	sp_buildtime = time()-buildtimestart
	
	#====================================================#

	#Solve IP
	status = optimize!(model)

	if termination_status(model) == MOI.OPTIMAL
		solvetime = solve_time(model)
		obj = objective_value(model)
		feasibleflag = 1
		#println("Total items picked = ", obj)

		if (debugmode == 1) & (debugprintstatements == 1)
			for w in sp.workstations, m in sp.orders, t in sp.times
				if value(v[m,w,t]) > 1e-4
					println("v[$m, $w, $t] = ", value(v[m,w,t]))
				end
			end
		end

	elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
		#println("Time out in the subproblem")
		solvetime = solve_time(model)
		obj = objective_value(model)
		#println("Total items picked = ", obj)
		feasibleflag = 0
	else
		println("Error in solving!")
		solvetime = solve_time(model)
		obj = 0
		#println("Total items picked = ", obj)
		feasibleflag = 0
	end

	#====================================================#

	return obj, solvetime, h, y, z, f, g, v, feasibleflag, sp_buildtime

end
