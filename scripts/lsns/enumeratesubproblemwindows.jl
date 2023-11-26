
using IterTools

function enumeratesubproblemwindows(partition, maxworkstationspersubproblem, subproblemtimelength)

	workstationgrouplist, timinglist = [], []

	allsubsets = []
	for sslength in 1:maxworkstationspersubproblem
		allsubsets = union(allsubsets, collect(subsets(partition.workstations,sslength)))
	end
	for wlist in allsubsets, tstart in -tstep:60:horizon+tstep-ceil((subproblemtimelength/length(wlist))/tstep)*tstep
		push!(workstationgrouplist, wlist)
		push!(timinglist, (tstart, tstart + ceil((subproblemtimelength/length(wlist))/tstep)*tstep))
	end

	#Initialize sets
	windows = []
	windowids = 1:length(workstationgrouplist)
	windowscontaining, windowsduring, windowidlookup = Dict(), Dict(), Dict()
	for w in partition.workstations, t in extendedtimes
		windowscontaining[w,t] = []
	end
	for t in extendedtimes
		windowsduring[t] = []
	end

	#Create window sets
	for sp in 1:length(workstationgrouplist)
		winobject = (workstations = workstationgrouplist[sp], times = max(0,timinglist[sp][1]):tstep:min(horizon,timinglist[sp][2]), extendedtimes = timinglist[sp][1]:tstep:timinglist[sp][2], tstart=timinglist[sp][1], tend=timinglist[sp][2])
		push!(windows, winobject)
		windowidlookup[winobject] = length(windows)

		for t in winobject.times
			push!(windowsduring[t], length(windows))
			for w in winobject.workstations
				push!(windowscontaining[w,t], length(windows))
			end
		end
	end

	return windows, windowsduring, windowidlookup, windowscontaining

end

