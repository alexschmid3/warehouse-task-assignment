
function updateglobalsolution(currpartition, currsol)

	for key in keys(currsol.y)
		globalsolution.y[key] = currsol.y[key]
	end
	for key in keys(currsol.h)
		globalsolution.h[key] = currsol.h[key]
	end
	for key in keys(currsol.v)
		globalsolution.v[key] = currsol.v[key]
	end
	for key in keys(currsol.ypath)
		globalsolution.ypath[key] = currsol.ypath[key]
	end

	for m in currpartition.orders
		remove!(globalsolution.unassignedorders, m)
	end
	for m in currsol.unassignedorders
		push!(globalsolution.unassignedorders, m)
	end

	for key in keys(currsol.podsworkedat)
		globalsolution.podsworkedat[key] = currsol.podsworkedat[key]
	end
	for key in keys(currsol.ordersopen)
		globalsolution.ordersopen[key] = currsol.ordersopen[key]
	end
	for key in keys(currsol.itempodpicklist)
		globalsolution.itempodpicklist[key] = currsol.itempodpicklist[key]
	end
	for key in keys(currsol.stationassign)
		globalsolution.stationassign[key] = currsol.stationassign[key]
	end
	for key in keys(currsol.podbusy)
		globalsolution.podbusy[key] = currsol.podbusy[key]
	end
	
end