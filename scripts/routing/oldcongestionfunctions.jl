

function getAplus(int1, int2, relevantintersections, leftstep, rightstep, upstep, downstep)

	moveright = intcoords[int2][1] - intcoords[int1][1] 
	movedown = intcoords[int2][2] - intcoords[int1][2] 

	A_plus = Dict()
	for i in relevantintersections[int1, int2]
		A_plus[i] = []
	end

	if moveright > 0
		for i in relevantintersections[int1, int2]
			if (intcoords[i][1] < intcoords[int2][1]) & (rightstep[i] != nothing)
				push!(A_plus[i], rightstep[i])
			end
		end
	elseif moveright < 0
		for i in relevantintersections[int1, int2]
			if (intcoords[i][1] > intcoords[int2][1]) & (leftstep[i] != nothing)
				push!(A_plus[i], leftstep[i])
			end
		end
	end

	if movedown > 0
		for i in relevantintersections[int1, int2]
			if (intcoords[i][2] < intcoords[int2][2]) & (downstep[i] != nothing)
				push!(A_plus[i], downstep[i])
			end
		end
	elseif movedown < 0
		for i in relevantintersections[int1, int2]
			if (intcoords[i][2] > intcoords[int2][2]) & (upstep[i] != nothing)
				push!(A_plus[i], upstep[i])
			end
		end
	end

	return A_plus

end

#-----------------------------------------------------------------------------------#

function explorepath(A_plus, originalint1, int1, int2, relevantintersections, traveltimeraw, currpath)

	for newint in intersect(relevantintersections[originalint1,int2], A_plus[int1])
		currentpath = deepcopy(currpath)
		push!(currentpath, newint)
		if newint != int2
			currentpath = explorepath(A_plus, originalint1, newint, int2, relevantintersections, traveltimeraw, currentpath)
		else
			push!(allpaths[originalint1,int2], currentpath)
			return currentpath
		end
	end

end

#-----------------------------------------------------------------------------------#

#int1, int2, relevantintersections, traveltimeraw, leftstep, rightstep, upstep, downstep, allpaths = int1, int2, relevantintersections, traveltimeraw, leftstep, rightstep, upstep, downstep, allpaths

function enumerateshortestpaths(int1, int2, relevantintersections, traveltimeraw, leftstep, rightstep, upstep, downstep)
	
	A_plus = getAplus(int1, int2, relevantintersections, leftstep, rightstep, upstep, downstep)
	
	currentpath = explorepath(A_plus, int1, int1, int2, relevantintersections, traveltimeraw, [int1])

	#for path in allpaths[int1,int2]
	#	pathwithtime = [(path[1],0)]
	#	for i in 2:length(path)
	#		push!(pathwithtime, (path[i], last(pathwithtime)[2] + traveltimeraw[path[i-1], path[i]]))
	#	end 
	#	push!(allpathswithtimes[int1,int2], pathwithtime)
	#end

end
