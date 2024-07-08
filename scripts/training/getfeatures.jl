
#-----------------------------------------------------------------------------------#

function getinstancefeatures_static()

	#timeused = podprocesstime * sum(sum(length(podspicked[w,t]) for t in times) for w in workstations) + itemprocesstime * sum(sum(length(itempodpicklist[w,t]) for t in times) for w in workstations)
 
	instance_features = Dict("lsnsiterationsbefore" => lsnsiterationsbefore,
							"dynamicmlpass" => dynamicmlpass, 
							"greedymethod" => orderprioritization, 
							"throughpututilization" => fullprobnewobj / (horizon * numworkstations / (podprocesstime + itemprocesstime)), 
							"timeutilization" => timeused / (horizon * numworkstations), 
							"podsperitem" => numpods * num_items_per_pod / num_unique_items , 
							"orders" => orders, 
							"itemson" => itemson, 
							"podswith" => podswith, 
							"pods" => pods, 
							"workstations" => workstations, 
							"tstep" => tstep, 
							"horizon" => horizon)

	return instance_features

end

#-----------------------------------------------------------------------------#

function getentityfeatures_static()
	 
	w_features = Dict("centrality" => Dict())
	t_features = Dict("horizonrelativetime" => Dict(), "timeofday" => Dict())
	p_features = Dict("avgdisttostations" => Dict(), "pctusefulitems" => Dict())
	m_features = Dict("ordersize" => Dict(), "oneitemflag" => Dict())
	i_features = Dict("iteminventory" => Dict(), "itemnumpods" => Dict())

	for w in workstations 
		w_features["centrality"][w] = sum(warehousedistance[l,w] for l in storagelocs) / length(storagelocs)
	end

	for t in extendedtimes
		t_features["horizonrelativetime"][t] = t / horizon
		t_features["timeofday"][t] = 0
	end

	for p in pods
		p_features["avgdisttostations"][p] = mean([manhattandist(podstorageloc[p],w) for w in workstations])
		p_features["pctusefulitems"][p] = length(podstartinventory[p]) / num_items_per_pod
	end

	for m in orders
		m_features["ordersize"][m] = length(itemson[m])
		m_features["oneitemflag"][m] = 2 - min(2, length(itemson[m]))
	end

	for i in items
		i_features["iteminventory"][i] = sum(inventory[i,p] for p in podswith[i])
		i_features["itemnumpods"][i] = length(podswith[i])
	end

	return w_features, t_features, p_features, m_features, i_features

end

#-----------------------------------------------------------------------------#

function getsynergyfeatures_static()

	mp_features = Dict("itemoverlap" => Dict(), "itemoverlappct" => Dict(), "overlapitemavginv" => Dict(), "overlapitemavgaltpods" => Dict())
	pw_features = Dict("distance" => Dict(), "xsteps" => Dict(), "ysteps" => Dict())
	pwt_features = Dict("existingcong" => Dict())
	wt_features = Dict("queuecong" => Dict(), "avglocalcong" => Dict(), "maxlocalcong" => Dict(), "avghyperlocalcong" => Dict(), "maxhyperlocalcong" => Dict())

	#Order-pod inventory features
	for m in orders, p in pods
		mp_features["itemoverlap"][m,p] = length(intersect(itemson[m], podstartinventory[p]))
		mp_features["itemoverlappct"][m,p] = mp_features["itemoverlap"][m,p] / length(itemson[m])
		if mp_features["itemoverlap"][m,p] > 1e-4
			mp_features["overlapitemavginv"][m,p] = mean([sum(inventory[i,p2] for p2 in podswith[i]) for i in intersect(itemson[m], podstartinventory[p])])
			mp_features["overlapitemavgaltpods"][m,p] = mean([length(podswith[i])-1 for i in intersect(itemson[m], podstartinventory[p])])
		else
			mp_features["overlapitemavginv"][m,p] = 0
			mp_features["overlapitemavgaltpods"][m,p] = 0
		end
	end

	#Pod-workstation-time congestion features
	#congestionat = Dict()
	#for l in reducedintersections, t in 0:congestiontstep:horizon
	#	if congestionarcs[l,t] != []
	#		congestionat[l,t] = sum(sum(congestioncontribution[a,l,t] * y_currsol[p,a] for p in arcpodset[a]) for a in congestionarcs[l,t])
	#	else 
	#		congestionat[l,t] = 0
	#	end
	#end
	relevantintersections = findrelevantintersection()
	for s in storagelocs, w in workstations, t in arclength[s,w]:tstep:horizon
		int1, int2 = maploctointersection[s], maploctointersection[w]
		allcongestionvals = []
		for l in intersect(relevantintersections[int1, int2], reducedintersections)
			localtraffic = congestionat[l,t] 
			push!(allcongestionvals, localtraffic / intersectionmaxpods[l])
		end
		for p in pods
			if podstorageloc[p] == s 
				pwt_features["existingcong"][p,w,t] = mean(allcongestionvals)
			end
		end
	end
	for s in storagelocs, w in workstations, t in 0:tstep:arclength[s,w]-tstep
		for p in pods
			if podstorageloc[p] == s 
				pwt_features["existingcong"][p,w,t] = 0
			end
		end
	end

	#Workstation-time congestion features
	kdtree = KDTree(transpose(intcoords_nn))
	mappingnum = reducedintersections[1] - 1
	for w in workstations, t in tstep:tstep:horizon
		#Queue congestion
		q = maploctointersection[w]
		wt_features["queuecong"][w,t] = mean([congestionat[q,t2] / intersectionmaxpods[q] for t2 in t-tstep:congestiontstep:t])
		
		#Local congestion
		closeints, dists = knn(kdtree, [intcoords[q][1], intcoords[q][2]], numlocalintersections, true)
		congestionnearworkstation = []
		for i in intersect([i2 + minimum(intersections) - 1 for i2 in closeints], reducedintersections)
			push!(congestionnearworkstation, mean([congestionat[i,t2] / intersectionmaxpods[i] for t2 in t-tstep:congestiontstep:t]))
		end
		wt_features["avglocalcong"][w,t] = mean(congestionnearworkstation)
		wt_features["maxlocalcong"][w,t] = maximum(congestionnearworkstation)

		#Hyper local congestion
		closestints, dists = knn(kdtree, [intcoords[q][1], intcoords[q][2]], numhyperlocalintersections, true)
		congestionnearworkstation = []
		for i in intersect([i2 + minimum(intersections) - 1 for i2 in closeints], reducedintersections)
			push!(congestionnearworkstation, mean([congestionat[i,t2] / intersectionmaxpods[i] for t2 in t-tstep:congestiontstep:t]))
		end
		wt_features["avghyperlocalcong"][w,t] = mean(congestionnearworkstation)
		wt_features["maxhyperlocalcong"][w,t] = maximum(congestionnearworkstation)
	end

	for w in workstations
		wt_features["queuecong"][w,0] = 0
		wt_features["avglocalcong"][w,0] = 0
		wt_features["maxlocalcong"][w,0] = 0
		wt_features["avghyperlocalcong"][w,0] = 0
		wt_features["maxhyperlocalcong"][w,0] = 0
	end

	#Pod-workstation features 
	for p in pods, w in workstations
		pw_features["distance"][p,w] = manhattandist(podstorageloc[p],w)
		pw_features["xsteps"][p,w] = abs(loccoords[podstorageloc[p],1] - loccoords[w,1]) 
		pw_features["ysteps"][p,w] = abs(loccoords[podstorageloc[p],2] - loccoords[w,2])
	end

	return mp_features, pw_features, pwt_features, wt_features

end

#-----------------------------------------------------------------------------#

function getcurrsolfeatures_static()
		
	currsol_features = Dict("mwt" => Dict(), "pwt" => Dict())

	for m in orders, w in workstations, t in tstep:tstep:horizon
		if v_currsol[m,w,t] > 0.001
			currsol_features["mwt"][m,w,t] = 1
		else
			currsol_features["mwt"][m,w,t] = 0
		end
	end

	for p in pods, w in workstations, t in 0:tstep:horizon
		if sum(y_currsol[p,a] for a in intersect(podarcset[p], union(A_plus[nodes[w,t]], A_minus[nodes[w,t]]))) > 0.001
			currsol_features["pwt"][p,w,t] = 1
		else
			currsol_features["pwt"][p,w,t] = 0
		end
	end

	openorders = Dict()
	for w in workstations, t in times
		openorders[w,t] = []
		for m in orders
			if v_currsol[m,w,t] > 0.001
				push!(openorders[w,t], m)
			end
		end
	end

	openpods = Dict()
	for w in workstations, t in times
		openpods[w,t] = []
		for p in pods
			if sum(y_currsol[p,a] for a in intersect(podarcset[p], union(A_plus[nodes[w,t]], A_minus[nodes[w,t]]))) > 0.001
				push!(openpods[w,t], p)
			end
		end
	end

	return currsol_features, openorders, openpods

end

#-----------------------------------------------------------------------------------#
