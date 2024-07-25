
function createroutingnodes()

	routenodelookup, routenodes, routenodes_end = Dict(), Dict(), []

	nindex = 1
	for t in 0:tstep_r:horizon_r, l in union(storagelocs, workstations)
		routenodes[l,t] = nindex
		routenodelookup[nindex] = (l,t)
		nindex += 1
	end
	numnodes_r = copy(nindex) - 1
	for t in union(-maxtraveltime:tstep_r:-tstep_r, horizon_r+tstep_r:tstep_r:horizon_r+maxtraveltime), l in union(storagelocs, workstations)
		routenodes[l,t] = nindex
		routenodelookup[nindex] = (l,t)
		if t == horizon_r+maxtraveltime
			push!(routenodes_end, nindex)
		end
		nindex += 1
	end
	extendednumnodes_r = copy(nindex) - 1

	return routenodelookup, routenodes, routenodes_end, numnodes_r, extendednumnodes_r

end

#-----------------------------------------------------------------------------------#

function findroutingprearcs()

	prearcs = []
	newarclength = Dict()
	#for l1 in setdiff(intersections, queueintersections), l2 in queueintersections
	for s in storagelocs, w in workstations
		l1, l2 = maploctointersection[s], maploctointersection[w] 
		dist = abs(intcoords[l1][1] - intcoords[l2][1]) + abs(intcoords[l1][2] - intcoords[l2][2]) 
		traveltime = dist / podspeed
		push!(prearcs, (s, w, traveltime, tstep_r * ceil(traveltime / tstep_r)))
		push!(prearcs, (w, s, traveltime, tstep_r * ceil(traveltime / tstep_r)))
		newarclength[w,s] = tstep_r * ceil(traveltime / tstep_r)
		newarclength[s,w] = tstep_r * ceil(traveltime / tstep_r)
	end

	return prearcs, newarclength

end

#-----------------------------------------------------------------------------------#

allpaths = Dict()

function pathsbetweenintersections()
	
	relevantintersections = getrelevantintersections() #GOOD
	traveltimeraw = getrawtraveltimes(podspeed)
	leftstep, rightstep, upstep, downstep = findpossiblesteps(horstreets, vertstreets, intlookup) #createintersectionarcs(horstreets, vertstreets, intlookup)

	for int1 in setdiff(intersections,queueintersections), int2 in queueintersections
		allpaths[int1,int2] = []
		enumerateshortestpaths(int1, int2, relevantintersections, traveltimeraw, leftstep, rightstep, upstep, downstep)
	end

end

#-----------------------------------------------------------------------------------#

function createroutingarcs()

	prearcs, newarclength = findroutingprearcs()
	pathsbetweenintersections()

	arcindex = 1
	routearclookup, routearclookup_path, routearcs, routearcs_path = Dict(), Dict(), Dict(), Dict()
	RA_plus, RA_minus, A_space = Dict(), Dict(), []
	for n in 1:extendednumnodes_r
		RA_plus[n] = []
		RA_minus[n] = []
	end

	#Physical arcs (incorporate paths)
	for pa in prearcs, t in -maxtraveltime:tstep_r:horizon_r+maxtraveltime
		l1, l2, tt = pa[1], pa[2], pa[4]
		if t + tt <= horizon_r+maxtraveltime
			n1, n2 = routenodes[l1, t], routenodes[l2, t + tt]
			routearcs[n1,n2] = []
			int1, int2 = maploctointersection[l1], maploctointersection[l2]
			if l1 != l2
				for rawpath in allpaths[min(int1, int2), max(int1, int2)]	
					if int1 == min(int1, int2)
						path = rawpath
					else
						path = reverse(rawpath)
					end
					push!(routearcs[n1, n2], arcindex)
					routearclookup[arcindex] = (n1, n2)
					routearclookup_path[arcindex] = (n1, n2, path)
					routearcs_path[n1, n2, path] = arcindex
					push!(RA_plus[n1], arcindex)
					push!(RA_minus[n2], arcindex)
					push!(A_space, arcindex)
					arcindex += 1
				end
			end
		end
	end

	#Stationary arcs
	for l in union(storagelocs, workstations), t in -maxtraveltime:tstep_r:horizon_r+maxtraveltime-tstep_r
		n1, n2 = routenodes[l, t], routenodes[l, t + tstep_r]
		routearcs[n1, n2] = [arcindex]
		path = [maploctointersection[l]]
		routearclookup[arcindex] = (n1, n2)
		routearclookup_path[arcindex] = (n1, n2, path)
		routearcs_path[n1, n2, path] = arcindex
		push!(RA_plus[n1], arcindex)
		push!(RA_minus[n2], arcindex)
		arcindex += 1
	end
	numarcs_r = arcindex - 1

	return routearcs, routearclookup, numarcs_r, RA_plus, RA_minus, newarclength, routearcs_path, routearclookup_path, A_space

end

#-----------------------------------------------------------------------------------#

function arcDesc(a)
	println(routenodelookup[routearclookup[a][1]], " --> ", routenodelookup[routearclookup[a][2]], " along : ", routearclookup_path[a][3])
end

#-----------------------------------------------------------------------------------#

function routingpodinfo()
	
	podstartnode_r = Dict()
	for p in pods_r
		podstartnode_r[p] = routenodes[podstorageloc[p], -maxtraveltime]
	end

	return podstartnode_r	

end

#-----------------------------------------------------------------------------------#

function podarcsets_routing(pods_r, newarclength)

	podarcset = Dict()
	A_minus_p, A_plus_p = Dict(), Dict()

	for p in pods_r
		podarcset[p] = []
	end
	for p in pods_r, n in 1:extendednumnodes_r
		A_plus_p[p,n] = []
		A_minus_p[p,n] = []
	end
	
	#Movement arcs
	for p in pods_r, w in workstations, t in -maxtraveltime:tstep_r:horizon_r+maxtraveltime-newarclength[podstorageloc[p], w]
		s = podstorageloc[p]
		n1_leave, n2_leave = routenodes[s, t], routenodes[w, t + newarclength[s, w]]
		n1_return, n2_return = routenodes[w, t], routenodes[s, t + newarclength[w, s]]

		for a_leave in routearcs[n1_leave, n2_leave]
			push!(podarcset[p], a_leave)
			push!(A_plus_p[p, n1_leave], a_leave)
			push!(A_minus_p[p, n2_leave], a_leave)
		end

		for a_return in routearcs[n1_return, n2_return]
			push!(podarcset[p], a_return)
			push!(A_plus_p[p, n1_return], a_return)
			push!(A_minus_p[p, n2_return], a_return)
		end
	end

	#Stationary arcs
	for p in pods_r, t in -maxtraveltime:tstep_r:horizon_r+maxtraveltime-tstep
		s = podstorageloc[p]
		for a in routearcs[routenodes[s, t], routenodes[s, t + tstep]]	
			push!(podarcset[p], a)
			#push!(arcpodset[a], p)
			push!(A_plus_p[p, routenodes[s, t]], a)
			push!(A_minus_p[p,routenodes[s, t + tstep]], a)
		end

		for w in workstations, a in routearcs[routenodes[w, t], routenodes[w, t + tstep]]
			push!(podarcset[p], a)
			#push!(arcpodset[a], p)
			push!(A_plus_p[p, routenodes[w, t]], a)
			push!(A_minus_p[p,routenodes[w, t + tstep]], a)
		end
	end
	
	return podarcset, A_minus_p, A_plus_p, podnodeset

end

#-----------------------------------------------------------------------------------#

function getpathcongestioncontributions()

	#Initialize
	congcontribarcs = Dict()
	for l in intersections, t in -maxtraveltime:tstep_cong:horizon_r+maxtraveltime
		congcontribarcs[l,t] = []
	end
	passingintersections, passingintersectiontimes = Dict(), Dict()
	for a in 1:numarcs_r
		passingintersections[a] = []
		passingintersectiontimes[a] = []
	end

	#Find the intersection-times for each arc
	for a in 1:numarcs_r
		n1, n2, path = routearclookup_path[a]
		t0, te = routenodelookup[n1][2], routenodelookup[n2][2]
		i0 = path[1]
		if (length(path) > 1) || !(intersect(path,queueintersections) == [])
			for i in path
				dist = abs(intcoords[i0][1] - intcoords[i][1]) + abs(intcoords[i0][2] - intcoords[i][2]) 
				traveltime = dist / podspeed
				tt_rdd = tstep_cong * ceil(traveltime / tstep_cong)
				push!(congcontribarcs[i, t0 + tt_rdd], a)
				push!(passingintersections[a], i)
				push!(passingintersectiontimes[a], (i,t0 + tt_rdd))
				
				#Add the congestion at the end of the arc
				if (i == last(path)) & (i in queueintersections)
					for t in t0 + tt_rdd + tstep_cong:tstep_cong:te-tstep_cong
						push!(congcontribarcs[i, t], a)
						push!(passingintersections[a], i)
						push!(passingintersectiontimes[a], (i,t))
					end
				end
			end
		end
	end

	return congcontribarcs, passingintersections, passingintersectiontimes

end

#-----------------------------------------------------------------------------------#




