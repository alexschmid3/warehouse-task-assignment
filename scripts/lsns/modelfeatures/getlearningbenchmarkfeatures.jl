
function getinstancefeatures_static()

	podsperitem = numpods * num_items_per_pod / num_unique_items
	instance_features = [0 1 1 podsperitem length(orders) length(allpods) length(workstations) tstep horizon]

	return instance_features

end

#-------------------------------------------------------------------------#

function preprocessentityfeatures_static()

	centrality, horizonrelativetime, avgdisttostations, ordersize, oneitemflag, iteminventory, itemnumpods = Dict(), Dict(), Dict(), Dict(), Dict(), Dict(), Dict()
	for w in workstations
		centrality[w] = sum(warehousedistance[l,w] for l in storagelocs) / length(storagelocs) 
	end
	for t in times
		horizonrelativetime[t] = t / horizon 
	end
	avgdisttostations = [mean([manhattandist(podstorageloc[p],w) for w in workstations]) for p in allpods]
	for m in orders 
		ordersize[m] = length(itemson[m])
		oneitemflag[m] = 2 - min(2, length(itemson[m])) 
	end
	for i in items
		iteminventory = sum(inventory[i,p] for p in podswith[i]) 
		itemnumpods[i] = length(podswith[i]) 
	end

	return centrality, horizonrelativetime, avgdisttostations, ordersize, oneitemflag, iteminventory, itemnumpods

end

#-------------------------------------------------------------------------#

function getentityfeatures_static(sp)

	sp_times_reg = max(sp.times[1], 0):tstep:min(last(sp.times), horizon)

	avg_centrality = mean([centrality[w] for w in sp.workstations])
	sum_centrality = sum([centrality[w] for w in sp.workstations])
	avg_horizonrelativetime = mean([horizonrelativetime[w] for w in sp_times_reg])
	sum_horizonrelativetime = sum([horizonrelativetime[w] for w in sp_times_reg])
	timeofday = 0
	avg_avgdisttostations = mean([avgdisttostations[w] for w in sp.pods])
	sum_avgdisttostations = sum([avgdisttostations[w] for w in sp.pods])
	avg_ordersize = mean([ordersize[w] for w in sp.orders])
	sum_ordersize = sum([ordersize[w] for w in sp.orders])
	pct_oneitemflag = mean([oneitemflag[w] for w in sp.orders])
	num_oneitemflag = sum([oneitemflag[w] for w in sp.orders])

	w_features = [avg_centrality sum_centrality]
	t_features = [avg_horizonrelativetime sum_horizonrelativetime timeofday]
	p_features = [avg_avgdisttostations sum_avgdisttostations]
	m_features = [avg_ordersize sum_ordersize pct_oneitemflag num_oneitemflag]

	return w_features, t_features, p_features, m_features

end

#-----------------------------------------------------------------------------#

function preprocesssynergyfeatures_static()

	itemoverlap, itemoverlappct, overlapitemavginv, overlapitemavgaltpods = Dict(), Dict(), Dict(), Dict()
	pwdistance, xsteps, ysteps = Dict(), Dict(), Dict()

	#Order-pod inventory features
	for m in orders, p in pods
		itemoverlap[m,p] = length(intersect(itemson[m], podstartinventory[p]))
		itemoverlappct[m,p] = itemoverlap[m,p] / length(itemson[m])
		if itemoverlap[m,p] > 1e-4
			overlapitemavginv[m,p] = mean([sum(inventory[i,p2] for p2 in podswith[i]) for i in intersect(itemson[m], podstartinventory[p])])
			overlapitemavgaltpods[m,p] = mean([length(podswith[i])-1 for i in intersect(itemson[m], podstartinventory[p])])
		else
			overlapitemavginv[m,p] = 0
			overlapitemavgaltpods[m,p] = 0
		end
	end

	#Pod-workstation features 
	for p in pods, w in workstations
		pwdistance[p,w] = manhattandist(podstorageloc[p],w)
		xsteps[p,w] = abs(loccoords[podstorageloc[p],1] - loccoords[w,1]) 
		ysteps[p,w] = abs(loccoords[podstorageloc[p],2] - loccoords[w,2])
	end

	return itemoverlap, itemoverlappct, overlapitemavginv, overlapitemavgaltpods, pwdistance, xsteps, ysteps

end

#-----------------------------------------------------------------------------#

function getsynergyfeatures_static(sp)

	avg_itemoverlap = sum(sum(itemoverlap[m,p] for p in sp.pods) for m in sp.orders) / sum(sum(1 for p in sp.pods) for m in sp.orders)
	sum_itemoverlap = sum(sum(itemoverlap[m,p] for p in sp.pods) for m in sp.orders)
	avg_itemoverlappct = sum(sum(itemoverlappct[m,p] for p in sp.pods) for m in sp.orders) / sum(sum(1 for p in sp.pods) for m in sp.orders)
	sum_itemoverlappct = sum(sum(itemoverlappct[m,p] for p in sp.pods) for m in sp.orders)
	avg_overlapitemavginv = sum(sum(overlapitemavginv[m,p] for p in sp.pods) for m in sp.orders) / sum(sum(1 for p in sp.pods) for m in sp.orders)
	sum_overlapitemavginv = sum(sum(overlapitemavginv[m,p] for p in sp.pods) for m in sp.orders)
	avg_overlapitemavgaltpods = sum(sum(overlapitemavgaltpods[m,p] for p in sp.pods) for m in sp.orders) / sum(sum(1 for p in sp.pods) for m in sp.orders)
	sum_overlapitemavgaltpods = sum(sum(overlapitemavgaltpods[m,p] for p in sp.pods) for m in sp.orders)

	mp_features = [avg_itemoverlap sum_itemoverlap avg_itemoverlappct sum_itemoverlappct avg_overlapitemavginv sum_overlapitemavginv avg_overlapitemavgaltpods sum_overlapitemavgaltpods]
	
	pw_total = sum(sum(1 for p in sp.pods) for w in sp.workstations)

	avg_distance = sum(sum(pwdistance[p,w] for p in sp.pods) for w in sp.workstations) / pw_total
	sum_distance = sum(sum(pwdistance[p,w] for p in sp.pods) for w in sp.workstations)
	avg_xsteps = sum(sum(xsteps[p,w] for p in sp.pods) for w in sp.workstations) / pw_total
	sum_xsteps = sum(sum(xsteps[p,w] for p in sp.pods) for w in sp.workstations)
	avg_ysteps = sum(sum(ysteps[p,w] for p in sp.pods) for w in sp.workstations) / pw_total
	sum_ysteps = sum(sum(ysteps[p,w] for p in sp.pods) for w in sp.workstations)
	
	pw_features = [avg_distance sum_distance avg_xsteps sum_xsteps avg_ysteps sum_ysteps]

	#Pod-workstation-time congestion features
	existingcong = Dict()
	relevantintersections = findrelevantintersection()
	for s in storagelocs, w in sp.workstations, t in max(sp.tstart,arclength[s,w]):tstep:min(sp.tend, horizon)
		int1, int2 = maploctointersection[s], maploctointersection[w]
		allcongestionvals = []
		for l in intersect(relevantintersections[int1, int2], reducedintersections)
			localtraffic = congestionat[l,t] 
			push!(allcongestionvals, localtraffic / intersectionmaxpods[l])
		end
		for p in sp.pods
			if podstorageloc[p] == s 
				existingcong[p,w,t] = mean(allcongestionvals)
			end
		end
	end
	for s in storagelocs, w in sp.workstations, t in max(0,sp.tstart):tstep:min(arclength[s,w]-tstep)
		for p in sp.pods
			if podstorageloc[p] == s 
				existingcong[p,w,t] = 0
			end
		end
	end

	avg_existingcong = sum(sum(sum(existingcong[p,w,t] for p in sp.pods) for t in max(0,sp.tstart):tstep:min(horizon,sp.tend)) for w in sp.workstations) / sum(sum(sum(1 for p in sp.pods) for t in max(0,sp.tstart):tstep:min(horizon,sp.tend)) for w in sp.workstations)
	sum_existingcong = sum(sum(sum(existingcong[p,w,t] for p in sp.pods) for t in max(0,sp.tstart):tstep:min(horizon,sp.tend)) for w in sp.workstations)

	pwt_features = [avg_existingcong sum_existingcong]

	#Workstation-time congestion features
	queuecong, avglocalcong, maxlocalcong, avghyperlocalcong, maxhyperlocalcong = Dict(), Dict(), Dict(), Dict(), Dict()
	kdtree = KDTree(transpose(intcoords_nn))
	mappingnum = reducedintersections[1] - 1
	for w in sp.workstations, t in max(sp.tstart, tstep):tstep:min(sp.tend, horizon)
		#Queue congestion
		q = maploctointersection[w]
		queuecong[w,t] = mean([congestionat[q,t2] / intersectionmaxpods[q] for t2 in t-tstep:congestiontstep:t])
		
		#Local congestion
		closeints, dists = knn(kdtree, [intcoords[q][1], intcoords[q][2]], numlocalintersections, true)
		congestionnearworkstation = []
		for i in intersect([i2 + minimum(intersections) - 1 for i2 in closeints], reducedintersections)
			push!(congestionnearworkstation, mean([congestionat[i,t2] / intersectionmaxpods[i] for t2 in t-tstep:congestiontstep:t]))
		end
		avglocalcong[w,t] = mean(congestionnearworkstation)
		maxlocalcong[w,t] = maximum(congestionnearworkstation)

		#Hyper local congestion
		closestints, dists = knn(kdtree, [intcoords[q][1], intcoords[q][2]], numhyperlocalintersections, true)
		congestionnearworkstation = []
		for i in intersect([i2 + minimum(intersections) - 1 for i2 in closeints], reducedintersections)
			push!(congestionnearworkstation, mean([congestionat[i,t2] / intersectionmaxpods[i] for t2 in t-tstep:congestiontstep:t]))
		end
		avghyperlocalcong[w,t] = mean(congestionnearworkstation)
		maxhyperlocalcong[w,t] = maximum(congestionnearworkstation)
	end

	for w in sp.workstations
		queuecong[w,0] = 0
		avglocalcong[w,0] = 0
		maxlocalcong[w,0] = 0
		avghyperlocalcong[w,0] = 0
		maxhyperlocalcong[w,0] = 0
	end

	avgtotal = sum(sum(1 for t in max(0,sp.tstart):tstep:min(horizon,sp.tend)) for w in sp.workstations)
	avg_queuecong = sum(sum(queuecong[w,t] for t in max(0,sp.tstart):tstep:min(horizon,sp.tend)) for w in sp.workstations) / avgtotal
	sum_queuecong = sum(sum(queuecong[w,t] for t in max(0,sp.tstart):tstep:min(horizon,sp.tend)) for w in sp.workstations)
	avg_avglocalcong = sum(sum(avglocalcong[w,t] for t in max(0,sp.tstart):tstep:min(horizon,sp.tend)) for w in sp.workstations) / avgtotal
	sum_avglocalcong = sum(sum(avglocalcong[w,t] for t in max(0,sp.tstart):tstep:min(horizon,sp.tend)) for w in sp.workstations)
	avg_maxlocalcong = sum(sum(maxlocalcong[w,t] for t in max(0,sp.tstart):tstep:min(horizon,sp.tend)) for w in sp.workstations) / avgtotal
	sum_maxlocalcong = sum(sum(maxlocalcong[w,t] for t in max(0,sp.tstart):tstep:min(horizon,sp.tend)) for w in sp.workstations)
	avg_avghyperlocalcong = sum(sum(avghyperlocalcong[w,t] for t in max(0,sp.tstart):tstep:min(horizon,sp.tend)) for w in sp.workstations) / avgtotal
	sum_avghyperlocalcong = sum(sum(avghyperlocalcong[w,t] for t in max(0,sp.tstart):tstep:min(horizon,sp.tend)) for w in sp.workstations)
	avg_maxhyperlocalcong = sum(sum(maxhyperlocalcong[w,t] for t in max(0,sp.tstart):tstep:min(horizon,sp.tend)) for w in sp.workstations) / avgtotal
	sum_maxhyperlocalcong = sum(sum(maxhyperlocalcong[w,t] for t in max(0,sp.tstart):tstep:min(horizon,sp.tend)) for w in sp.workstations)

	wt_features = [avg_queuecong sum_queuecong avg_avglocalcong sum_avglocalcong avg_maxlocalcong sum_maxlocalcong avg_avghyperlocalcong sum_avghyperlocalcong avg_maxhyperlocalcong sum_maxhyperlocalcong]

	return mp_features, pw_features, pwt_features, wt_features

end

#-------------------------------------------------------------------------#

function capacityslackfeatures(sp, currsol)

	sp_times_reg = max(sp.times[1], 0):tstep:min(last(sp.times), horizon)

	#Slack in capacity constraints
	stationitemspicked = sum(sum(sum(sum(sum(h_currsol[m,i,p,w,t] - h_currsol[m,i,p,w,t-tstep] for p in sp.podswith[i]) for i in sp.itemson[m]) for m in orders if sp.itemson[m] != []) for t in setdiff(sp.times, sp.times[1])) for w in sp.workstations) + sum(sum(sum(sum(h_currsol[m,i,p,w,sp_times[1]] for p in sp.podswith[i]) for i in sp.itemson[m]) for m in sp.orders if sp.itemson[m] != []) for w in sp.workstations)
	stationpodspicked = sum(sum(length(currsol.podsworkedat[w,t]) for w in sp.workstations) for t in sp_times_reg)
	stationmissingthroughput = sum(sum(tstep for t in sp_times_reg) for w in sp_workstations) - stationitemspicked*itemprocesstime - stationpodspicked*podprocesstime
	stationidlepct = stationmissingthroughput / sum(sum(tstep for t in sp_times_reg) for w in sp_workstations)
	orderslottimessopen = sum(sum(C[w] for t in sp_times_reg) for w in sp_workstations) - sum(sum(sum(v_currsol[m,w,t] for m in orders) for t in sp_times_reg) for w in sp_workstations)  

	capacityslackparms = [stationmissingthroughput stationidlepct orderslottimessopen]

	return capacityslackparms, stationitemspicked

end

#-------------------------------------------------------------------------#

function compatibilityfeatures(sp_workstations, sp_times, sp_orders)

	sp_times_reg = max(sp_times[1], 0):tstep:min(last(sp_times), horizon)

	#"Synergy" with other orders that could fill in gaps
	if (unassignedorders == []) || (setdiff(sp_orders, unassignedorders) == []) 
		podcompatibilitywithunassignedorders = 0
	else
		podcompatibilitywithunassignedorders = 0
		for m1 in setdiff(sp_orders, unassignedorders), m2 in intersect(sp_orders, unassignedorders)
			podcompatibilitywithunassignedorders += orderordercompat[m1, m2] 
		end
	end
	easyoneitemorders = 0
	for t in sp_times_reg, w in sp_workstations
		if (podsworkedat[w,t] != []) & (unassignedorders != [])
			for m in unassignedorders, p in podsworkedat[w,t]
				easyoneitemorders += orderpodcompat[m,p] 
			end
		end
	end	
	compatibilityparms = [podcompatibilitywithunassignedorders easyoneitemorders]

	return compatibilityparms

end

#-------------------------------------------------------------------------#

function stationstatusfeatures(sp_workstations, sp_times, sp_tstart, sp_tend, sp_orders, sp_itemson, subprobArcSet, y_known, sp_workedorders, stationitemspicked)

	sp_times_reg = max(sp_times[1], 0):tstep:min(last(sp_times), horizon)

	#Station current stats 
	totalworkedorders = length(sp_workedorders)
	if sp_workedorders != []
		stationavgordersize = sum(length(itemson[m]) for m in sp_workedorders) / totalworkedorders
		stationavgorderinventory = 0 #sum(sum(sum(inventory[i,p] for p in podswith[i]) for i in sp_itemson[m]) for m in sp_workedorders if sp_itemson[m] != []) / sum(length(sp_itemson[m]) for m in sp_workedorders if sp_itemson[m] != []) 								
	else
		stationavgordersize = 0
		stationavgorderinventory = 0
	end
	stationpodsvisited = 0 # sum(sum(sum(sum(y_currsol[p,a] for a in setdiff(intersect(A_plus_p[p,nodes[w,t]], subprobArcSet), A_queues)) for t in setdiff(sp_times_reg,last(sp_times_reg)) if setdiff(intersect(A_plus_p[p,nodes[w,t]], subprobArcSet), A_queues) != []) for w in sp_workstations) for p in pods) + sum(sum(sum(y_known[p,nodes[w,t]] for w in sp_workstations) for t in sp_times_reg) for p in pods)
	avgpicksperpod = 0 #stationitemspicked / stationpodsvisited
	#podavgdistancetraveled = sum(sum(sum(abs.(loccoords[w,:] - loccoords[podstorageloc[p],:])) * (sum(sum(y[p,a] for a in setdiff(intersect(A_plus_p[p,nodes[w,t]], subprobArcSet), A_queues)) for t in setdiff(sp_times, last(sp_times))) + y_known[p,nodes[w,last(sp_times)]]) for w in sp_workstations) for p in pods)	/ stationpodsvisited
	stationcurrentstatparms = [totalworkedorders stationavgordersize stationavgorderinventory stationpodsvisited avgpicksperpod]

	return stationcurrentstatparms

end

#-------------------------------------------------------------------------#

function algorithmcontolfeatures(subproblem_id, lsnsiter)

	#Algorithm control
	lsnsiterationssofar = lsnsiter
	recentlyoptimized_flag = windowrecentlyreoptimized[subproblem_id]
	totalsubproblemreoptimizations = windowreoptimizedcount[subproblem_id]
	mostrecentimprovement = mostrecentwindowimprovement[subproblem_id]

	algcontrolparms = [lsnsiterationssofar recentlyoptimized_flag totalsubproblemreoptimizations mostrecentimprovement]
	
	return algcontrolparms

end 

#-------------------------------------------------------------------------#

function calculateadditionalfeatures(subproblem_id, lsnsiter, sp_workstations, sp_times, sp_tstart, sp_tend, sp_orders, sp_itemson, ambientcongestion, subprobArcSet, y_known, assignedorders)

	capacityslackparms, stationitemspicked = capacityslackfeatures(sp_workstations, sp_times, sp_orders, sp_itemson, subprobArcSet)
	compatibilityparms = compatibilityfeatures(sp_workstations, sp_times, sp_orders)
	stationcurrentstatparms = stationstatusfeatures(sp_workstations, sp_times, sp_tstart, sp_tend, sp_orders, sp_itemson, subprobArcSet, y_known, assignedorders, stationitemspicked)
	algcontrolparms = algorithmcontolfeatures(subproblem_id, lsnsiter)

	features = hcat(capacityslackparms, compatibilityparms, stationcurrentstatparms, algcontrolparms)

	return features

end

#-------------------------------------------------------------------------#

function getsubproblemfeatures(sp, lsnsiter)

    id_list = [run_id instance_id warehouse_id sp.winid]
	static_features = getinstancefeatures_static()
	w_features, t_features, p_features, m_features = getentityfeatures_static(sp)
	mp_features, pw_features, pwt_features, wt_features = getsynergyfeatures_static(sp)
	additional_features = calculateadditionalfeatures(subproblem_id, lsnsiter, sp_workstations, sp_times, sp_tstart, sp_tend, sp_orders, sp_itemson, ambientcongestion, subprobArcSet, y_known, assignedorders)

	return hcat(id_list, static_features, w_features, t_features, p_features, m_features, mp_features, pw_features, pwt_features, wt_features, additional_features)

end