
function getorderpodlists()
	
	podsfor = Dict()

	for m in orders
		podsfor[m] = []
		for i in itemson[m]
			podsfor[m] = vcat(podsfor[m], podswith[i])
		end
	end

	return podsfor

end

#-------------------------------------------------------------------------#

function calculateordercompatibility()
	
	podsfor = getorderpodlists()
	compat = Dict()

	for i1 in 1:length(orders), i2 in i1+1:length(orders)
		m1, m2 = orders[i1], orders[i2]
		compat[m1,m2] = length(intersect(podsfor[m1], podsfor[m2]))
		compat[m2,m1] = compat[m1,m2]
	end

	return compat

end

#-------------------------------------------------------------------------#

function calculateorderpodcompatibility()

	compat = Dict()

	for m in orders, p in pods
		compat[m,p] = length(intersect(itemson[m], podstartinventory[p]))
	end

	return compat

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

#-------------------------------------------------------------------------#

function preprocesslearningbenchmarkfeatures()

    orderordercompat = calculateordercompatibility()
    orderpodcompat = calculateorderpodcompatibility()
    centrality, horizonrelativetime, avgdisttostations, ordersize, oneitemflag, iteminventory, itemnumpods = preprocessentityfeatures_static()
    itemoverlap, itemoverlappct, overlapitemavginv, overlapitemavgaltpods, pwdistance, xsteps, ysteps = preprocesssynergyfeatures_static()

    lbfeatures = (orderordercompat=orderordercompat, orderpodcompat=orderpodcompat, centrality=centrality, 
                horizonrelativetime=horizonrelativetime, avgdisttostations=avgdisttostations, ordersize=ordersize, 
                oneitemflag=oneitemflag, iteminventory=iteminventory, itemnumpods=itemnumpods, itemoverlap=itemoverlap, 
                itemoverlappct=itemoverlappct, overlapitemavginv=overlapitemavginv, overlapitemavgaltpods=overlapitemavgaltpods, 
                pwdistance=pwdistance, xsteps=xsteps, ysteps=ysteps)

    return lbfeatures

end