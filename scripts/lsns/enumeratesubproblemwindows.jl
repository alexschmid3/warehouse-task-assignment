
using IterTools

function enumeratesubproblemwindows(currpartition, maxworkstationspersubproblem, subproblemtimelength)

	workstationgrouplist, timinglist = [], []

	allsubsets = []
	#targetnumworkstations
	if methodparamsfilename == "data/extensions/spsize/test_run_parameters.csv"
		allsubsets = union(allsubsets, collect(subsets(currpartition.workstations,targetnumworkstations)))
		println("Station subsets = ", allsubsets)
	else
		for sslength in 1:maxworkstationspersubproblem
			allsubsets = union(allsubsets, collect(subsets(currpartition.workstations,sslength)))
		end
		println("Station subsets = ", allsubsets)
	end
	if warehouseparamsfilename == "data/extensions/timedisc/warehouse_sizes_and_capacities.csv"
		for wlist in allsubsets, tstart in -tstep:tstep:horizon+tstep-ceil((subproblemtimelength/length(wlist))/tstep)*tstep
			push!(workstationgrouplist, wlist)
			push!(timinglist, (tstart, tstart + ceil((subproblemtimelength/length(wlist))/tstep)*tstep))
		end
	else
		for wlist in allsubsets, tstart in -tstep:60:horizon+tstep-ceil((subproblemtimelength/length(wlist))/tstep)*tstep
			push!(workstationgrouplist, wlist)
			push!(timinglist, (tstart, tstart + ceil((subproblemtimelength/length(wlist))/tstep)*tstep))
		end
	end

	#Initialize sets
	windows = []
	windowids = 1:length(workstationgrouplist)
	windowscontaining, windowsduring, windowidlookup = Dict(), Dict(), Dict()
	for w in currpartition.workstations, t in extendedtimes
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

