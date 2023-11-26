
function getreducedintersections2(partition)
	
	max_y_coord = maximum([intcoords[i][2] for i in partition.intersections if !(i in queueintersections)])
	min_y_coord = minimum([intcoords[i][2] for i in partition.intersections if !(i in queueintersections)])
	min_x_coord = minimum([intcoords[i][1] for i in partition.intersections if !(i in queueintersections)])
	max_x_coord = maximum([intcoords[i][1] for i in partition.intersections if !(i in queueintersections)])
	x_interval = 8
	y_interval = 6

	reducedintersections = []
	for i in partition.intersections
		intx, inty = intcoords[i]
		if inty >= max_y_coord
			push!(reducedintersections, i)
		elseif mod((intx - min_x_coord)/ x_interval, 2) == mod((inty - min_y_coord)/ y_interval, 2)  
			push!(reducedintersections, i)
		elseif i in queueintersections
			push!(reducedintersections, i)
		end
	end

	return reducedintersections

end

#-----------------------------------------------------------------------------------------------------#

function findrelevantintersection2(partition)

	relevantintersections = Dict() 

	for int1 in setdiff(partition.intersections,queueintersections), int2 in intersect(partition.intersections,queueintersections)
		relevantintersections[int1, int2] = []
		for i in union(setdiff(partition.intersections, queueintersections), int2)
			if (intcoords[i][1] >= min(intcoords[int1][1], intcoords[int2][1])) & (intcoords[i][1] <= max(intcoords[int1][1], intcoords[int2][1])) & (intcoords[i][2] >= min(intcoords[int1][2], intcoords[int2][2])) & (intcoords[i][2] <= max(intcoords[int1][2], intcoords[int2][2]))
				push!(relevantintersections[int1, int2], i)
			end
		end
	end
	
	return relevantintersections

end

#-----------------------------------------------------------------------------------------------------#

function intersectiontraveltimes2(partition, podspeed)

	traveltimeraw = Dict()

	for int1 in partition.intersections, int2 in partition.intersections
		manhattandist2 = abs(intcoords[int1][1] - intcoords[int2][1]) +  abs(intcoords[int1][2] - intcoords[int2][2])
		arclengthraw = manhattandist2 / podspeed
		traveltimeraw[int1, int2] = arclengthraw
	end

	return traveltimeraw

end

#-----------------------------------------------------------------------------------------------------#

function createintersectionarcs2(horstreets, vertstreets, intlookup)
	
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

#-----------------------------------------------------------------------------------------------------#

function directcontributioncalculation2(partition)
	
	relevantintersections = findrelevantintersection2(partition) #GOOD
	traveltimeraw = intersectiontraveltimes2(partition, podspeed)
	leftstep, rightstep, upstep, downstep = createintersectionarcs2(horstreets, vertstreets, intlookup)

	trafficcontribution = Dict()
	for int1 in setdiff(partition.intersections, queueintersections), int2 in intersect(partition.intersections, queueintersections)
		for int3 in setdiff(relevantintersections[int1,int2], int2)
			interimint2 = intlookup[intcoords[int2][1], intcoords[int2][2]-3]
			trafficcontribution[int1,int2,int3] = stepbasedcontributioncalc(int1,interimint2,int3)
		end
		trafficcontribution[int1,int2,int2] = 1.0
	end

	return trafficcontribution, traveltimeraw

end

#-----------------------------------------------------------------------------------------------------#

function directcongestioncalculation2(partition, trafficcontribution, reducedintersections, traveltimeraw)
	
	relevantintersections = findrelevantintersection()

	congestioncontribution = Dict()
	congestionarcs = Dict()
	for l in partition.intersections, t in 0:congestiontstep:horizon
		congestionarcs[l,t] = []
	end

	#Add movement arcs to congestion
	for s in partition.storagelocs, ws in partition.workstations
		int1, int2 = maploctointersection[s], maploctointersection[ws]
		for t in 0:tstep:horizon-arclength[s,ws]
			a1 = arcs[nodes[s, t], nodes[ws, t + arclength[s,ws]]]
			a2 = arcs[nodes[ws, t], nodes[s, t + arclength[s,ws]]]
			for int3 in intersect(relevantintersections[int1, int2], reducedintersections)
				tt_in = convert(Int64, floor(traveltimeraw[int1, int3]/congestiontstep) * congestiontstep)
				tt_out = convert(Int64, floor(traveltimeraw[int3, int2]/congestiontstep) * congestiontstep)
				congestioncontribution[a1, int3, t+tt_in] = trafficcontribution[int1, int2, int3]
				congestioncontribution[a2, int3, t+tt_out] = trafficcontribution[int1, int2, int3]
				push!(congestionarcs[int3, t+tt_in], a1)
				push!(congestionarcs[int3, t+tt_out], a2)
			end
		end
	end

	#Add queue congestion for stationary arcs
	for ws in partition.workstations, t in 0:tstep:horizon-tstep
		a = arcs[nodes[ws, t], nodes[ws, t+tstep]]
		int3 = maploctointersection[ws]
		for t2 in t:congestiontstep:t+tstep-congestiontstep
			congestioncontribution[a, int3, t2] = 1
			push!(congestionarcs[int3, t2], a)
		end
	end

	return congestioncontribution, congestionarcs

end

#-----------------------------------------------------------------------------------------------------#

#include("scripts/congestionmodel_decomp.jl")

#-----------------------------------------------------------------------------------------------------#

function configurepartitioncongestion(partition)

	reducedintersections = getreducedintersections2(partition)
	trafficcontribution, traveltimeraw = directcontributioncalculation2(partition)
	congestioncontribution, congestionarcs = directcongestioncalculation2(partition, trafficcontribution, reducedintersections, traveltimeraw)

	partitionfloor = (congestioncontribution=congestioncontribution, congestionarcs=congestionarcs, trafficcontribution=trafficcontribution, 
					traveltimeraw=traveltimeraw, reducedintersections=reducedintersections)

	return partitionfloor 

end
