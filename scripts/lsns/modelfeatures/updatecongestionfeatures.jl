
#-----------------------------------------------------------------------------#

function updatecongestionfeatures()
	
	#Congestion pre-processing
	allcong = sum(currcong[p] for p in pods)
	maxcong = zeros(length(intersections), length(0:congestiontstep:horizon))
	for i in intersections, t in 1:length(0:congestiontstep:horizon)
		maxcong[maps.mapintersectiontorow[i], t] = intersectionmaxpods[i]
	end
	normalizedcong = allcong ./ maxcong

	#Workstation-time congestion features
	kdtree = KDTree(transpose(intcoords_nn))
	for w in workstations, t in 0:tstep:horizon
		#Queue congestion
		q = maps.mapintersectiontorow[maploctointersection[w]]
		features.wt["queuecong"][w-last(storagelocs),convert(Int,t/tstep+1)] = mean([normalizedcong[q,maps.maptimetocolumn[t2]] for t2 in max(0,t-tstep):congestiontstep:t])
	
		#Local congestion
		closeints, dists = knn(kdtree, [intcoords[q][1], intcoords[q][2]], numlocalintersections, true)
		congestionnearworkstation = [sum(normalizedcong[maps.mapintersectiontorow[i+minimum(intersections)-1],maps.maptimetocolumn[t2]] for t2 in max(0,t-tstep):congestiontstep:t) for i in closeints]
		features.wt["avglocalcong"][w-last(storagelocs),convert(Int,t/tstep+1)] = mean(congestionnearworkstation)
		features.wt["maxlocalcong"][w-last(storagelocs),convert(Int,t/tstep+1)] = maximum(congestionnearworkstation)

		#Hyper local congestion
		closestints, dists = knn(kdtree, [intcoords[q][1], intcoords[q][2]], numhyperlocalintersections, true)
		congestionnearestworkstation = [sum(normalizedcong[maps.mapintersectiontorow[i+minimum(intersections)-1],maps.maptimetocolumn[t2]] for t2 in max(0,t-tstep):congestiontstep:t) for i in closestints]
		features.wt["avghyperlocalcong"][w-last(storagelocs),convert(Int,t/tstep+1)] = mean(congestionnearestworkstation)
		features.wt["maxhyperlocalcong"][w-last(storagelocs),convert(Int,t/tstep+1)] = maximum(congestionnearestworkstation)
	end

end

#-----------------------------------------------------------------------------#