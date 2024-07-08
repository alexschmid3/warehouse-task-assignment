
function findimpactordersandpods(sp_workstations, sp_orders, sp_pods, sp_itemson, sp_times, sp_podswith, h_currsol, h_sp)

	impact_orders, impact_pods = [], []
	lasttime = min(horizon, maximum(sp_times))

	for m in sp_orders
		try
			if sp_itemson[m] != []
				newitemspicked = sum(sum(sum(value(h_sp[m,i,p,w,lasttime]) for p in sp_podswith[i]) for i in sp_itemson[m] if sp_podswith[i] != []) for w in sp_workstations)
				if newitemspicked > 0.001
					push!(impact_orders, m)
				end
			end
		catch
			1+1
		end
	end
	for p in sp_pods
		newpodpicks = 0
		for w in sp_workstations, i in allitems 
			if p in podswith[i]
				for m in sp_orders 
					if i in sp_itemson[m]
						newpodpicks += value(h_sp[m,i,p,w,lasttime])
					end
				end
			end
		end
		if newpodpicks > 0.001
			push!(impact_pods, p)
		end
	end

	return unique!(impact_orders), unique!(impact_pods)

end
