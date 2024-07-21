
function writeordergraph(ordergraphreportingfilename)

	mystationlist, mytimelist, myorderlist, mypodlist, myitemlist = [], [], [], [], []

	for s in 1:numpartitions
		currpartition, currsol = partitioninfo[s], partitionsolution[s]
		for w in currpartition.workstations, t in times, (m,i,p) in currsol.itempodpicklist[w,t]
			push!(mystationlist, w)
			push!(mytimelist, t)
			push!(myorderlist, m)
			push!(mypodlist, p)
			push!(myitemlist, i)
		end
	end

	df = DataFrame(station=mystationlist, time=mytimelist, order=myorderlist, pod=mypodlist, item=myitemlist)

	CSV.write(ordergraphreportingfilename, df)

end
