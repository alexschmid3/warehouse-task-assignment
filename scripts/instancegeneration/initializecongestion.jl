
using SparseArrays

#---------------------------------------------------------------------------------------#

function emptycongestion()

	currentcongestion = Dict()

	for p in pods
		currentcongestion[p] = spzeros(length(intersections), length(0:congestiontstep:horizon))
	end

	return currentcongestion

end

#---------------------------------------------------------------------------------------#

function createcongestionmaps()

	mapcolumntotime, maprowtointersection, maptimetocolumn, mapintersectiontorow = Dict(), Dict(), Dict(), Dict()
	rowindex, colindex = 1, 1
	for i in intersections
		mapintersectiontorow[i] = rowindex
		maprowtointersection[rowindex] = i
		rowindex += 1
	end
	for t in 0:congestiontstep:horizon
		maptimetocolumn[t] = colindex
		mapcolumntotime[colindex] = t
		colindex += 1
	end

	maps = (mapcolumntotime=mapcolumntotime, maprowtointersection=maprowtointersection, maptimetocolumn=maptimetocolumn, mapintersectiontorow=mapintersectiontorow)

	return maps

end

#---------------------------------------------------------------------------------------#

function getrawtraveltimes(podspeed)

	traveltimeraw = Dict()

	for int1 in intersections, int2 in intersections
		manhattandist2 = abs(intcoords[int1][1] - intcoords[int2][1]) +  abs(intcoords[int1][2] - intcoords[int2][2])
		arclengthraw = manhattandist2 / podspeed
		traveltimeraw[int1, int2] = arclengthraw
	end

	return traveltimeraw

end

#---------------------------------------------------------------------------------------#

function findpossiblesteps(horstreets, vertstreets, intlookup)
	
	#intarcs, intA_plus, intA_minus = [], Dict(), Dict()
	leftstep, rightstep, upstep, downstep = Dict(), Dict(), Dict(), Dict()

	for i in intersections
		leftstep[i] = nothing
		rightstep[i] = nothing
		upstep[i] = nothing
		downstep[i] = nothing
	end

	#Right steps
	for hs in 1:size(horstreets)[1], vs1 in 1:size(vertstreets)[1]-1, 
		vs2 = vs1+1
		inty = (horstreets[hs,1] + horstreets[hs,2]) / 2
		intx1 = (vertstreets[vs1,1] + vertstreets[vs1,2]) / 2
		intx2 = (vertstreets[vs2,1] + vertstreets[vs2,2]) / 2
		rightstep[intlookup[intx1, inty]] = intlookup[intx2, inty]
	end

	#Left steps
	for hs in 1:size(horstreets)[1], vs1 in 2:size(vertstreets)[1]
		vs2 = vs1-1
		inty = (horstreets[hs,1] + horstreets[hs,2]) / 2
		intx1 = (vertstreets[vs1,1] + vertstreets[vs1,2]) / 2
		intx2 = (vertstreets[vs2,1] + vertstreets[vs2,2]) / 2
		leftstep[intlookup[intx1, inty]] = intlookup[intx2, inty]
	end

	#Up steps
	for vs in 1:size(vertstreets)[1], hs1 in 2:size(horstreets)[1], 
		hs2 = hs1-1
		inty1 = (horstreets[hs1,1] + horstreets[hs1,2]) / 2
		inty2 = (horstreets[hs2,1] + horstreets[hs2,2]) / 2
		intx = (vertstreets[vs,1] + vertstreets[vs,2]) / 2
		upstep[intlookup[intx, inty1]] = intlookup[intx, inty2]
	end

	#Down steps
	for vs in 1:size(vertstreets)[1], hs1 in 1:size(horstreets)[1]-1, 
		hs2 = hs1+1
		inty1 = (horstreets[hs1,1] + horstreets[hs1,2]) / 2
		inty2 = (horstreets[hs2,1] + horstreets[hs2,2]) / 2
		intx = (vertstreets[vs,1] + vertstreets[vs,2]) / 2
		downstep[intlookup[intx, inty1]] = intlookup[intx, inty2]
	end

	#All steps for queue intersections
	for ws in workstations
		intw = maploctointersection[ws]
		ints = intlookup[intcoords[intw][1], intcoords[intw][2]-3]
		upstep[intw] = ints
		downstep[ints] = intw
	end

	#return intarcs, intA_plus, intA_minus
	return leftstep, rightstep, upstep, downstep

end

#---------------------------------------------------------------------------------------#

function getrelevantintersections()

	relevantintersections = Dict() 

	for int1 in setdiff(intersections,queueintersections), int2 in intersect(intersections,queueintersections)
		relevantintersections[int1, int2] = []
		for i in union(setdiff(intersections, queueintersections), int2)
			if (intcoords[i][1] >= min(intcoords[int1][1], intcoords[int2][1])) & (intcoords[i][1] <= max(intcoords[int1][1], intcoords[int2][1])) & (intcoords[i][2] >= min(intcoords[int1][2], intcoords[int2][2])) & (intcoords[i][2] <= max(intcoords[int1][2], intcoords[int2][2]))
				push!(relevantintersections[int1, int2], i)
			end
		end
	end
	
	return relevantintersections

end

#---------------------------------------------------------------------------------------#

function stepbasedcontributioncalc_refactor(int1,int2,int3)

	moveright = convert(Int, abs( intcoords[int2][1] - intcoords[int1][1] ) / 8 )
	movedown = convert(Int, abs( intcoords[int2][2] - intcoords[int1][2] ) / 6 )

	moveright_before = convert(Int, abs( intcoords[int3][1] - intcoords[int1][1] ) / 8 )
	movedown_before = convert(Int, abs( intcoords[int3][2] - intcoords[int1][2] ) / 6 )

	moveright_after = convert(Int, abs( intcoords[int2][1] - intcoords[int3][1] ) / 8 ) 
	movedown_after = convert(Int, abs( intcoords[int2][2] - intcoords[int3][2] ) / 6 )

	waystogettoint3 = factorial(big(moveright_before + movedown_before)) / (factorial(big(moveright_before)) * factorial(big(movedown_before)))
	waystogetfromint3 = factorial(big(moveright_after + movedown_after)) / (factorial(big(moveright_after)) * factorial(big(movedown_after)))

	pathsthroughint3 = waystogettoint3 * waystogetfromint3
	totalpaths = factorial(big(moveright + movedown)) / (factorial(big(moveright)) * factorial(big(movedown)))

	return pathsthroughint3 / totalpaths

end

#---------------------------------------------------------------------------------------#

function trafficcontributioncalc()
	
	relevantintersections = getrelevantintersections()
	traveltimeraw = getrawtraveltimes(podspeed)
	leftstep, rightstep, upstep, downstep = findpossiblesteps(horstreets, vertstreets, intlookup)

	trafficcontribution = Dict()
	for int1 in setdiff(intersections, queueintersections), int2 in intersect(intersections, queueintersections)
		for int3 in setdiff(relevantintersections[int1,int2], int2)
			interimint2 = intlookup[intcoords[int2][1], intcoords[int2][2]-3]
			trafficcontribution[int1,int2,int3] = stepbasedcontributioncalc_refactor(int1,interimint2,int3)
		end
		trafficcontribution[int1,int2,int2] = 1.0
	end

	return trafficcontribution, traveltimeraw, relevantintersections

end

#---------------------------------------------------------------------------------------#

function createcongestionsignatures(maps)

	trafficcontribution, traveltimeraw, relevantintersections = trafficcontributioncalc()

	congestionsignature = Dict()
	for a in 1:numarcs
		congestionsignature[a] = spzeros(length(intersections), length(0:congestiontstep:horizon))
	end

	for s in storagelocs, ws in workstations
		int1, int2 = maploctointersection[s], maploctointersection[ws]
		for t in 0:tstep:horizon-arclength[s,ws]
			a1 = arcs[nodes[s, t], nodes[ws, t + arclength[s,ws]]]
			a2 = arcs[nodes[ws, t], nodes[s, t + arclength[s,ws]]]

			for int3 in intersect(relevantintersections[int1, int2], intersections)
				tt_in = convert(Int64, floor(traveltimeraw[int1, int3]/congestiontstep) * congestiontstep)
				tt_out = convert(Int64, floor(traveltimeraw[int3, int2]/congestiontstep) * congestiontstep)

				int3_index = maps.mapintersectiontorow[int3]
				t1_index = maps.maptimetocolumn[t+tt_in]
				t2_index = maps.maptimetocolumn[t+tt_out]

				congestionsignature[a1][int3_index, t1_index] = trafficcontribution[int1, int2, int3]
				if !(int3 in queueintersections)
					congestionsignature[a2][int3_index, t2_index] = trafficcontribution[int1, int2, int3]
				end
			end
		end
	end

	return congestionsignature

end

#---------------------------------------------------------------------------------------#

function getintersectionmaxpods()

	intersectionmaxpods = Dict()
	for i in intersections
		if i in queueintersections
			intersectionmaxpods[i] = intersectioncapacity * 2
		else
			intersectionmaxpods[i] = intersectioncapacity
		end
	end

	return intersectionmaxpods

end

#---------------------------------------------------------------------------------------#

function initializecongestion()

	currentcongestion = emptycongestion()
	maps = createcongestionmaps()
	congestionsignature = createcongestionsignatures(maps)
	intersectionmaxpods = getintersectionmaxpods()

	return currentcongestion, maps, congestionsignature, intersectionmaxpods

end
