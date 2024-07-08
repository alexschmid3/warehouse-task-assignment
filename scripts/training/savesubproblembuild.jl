
function savesubproblembuild(subproblem_id, sp_workstations, sp_times, impact_orders, impact_pods, new_obj, old_obj)

	#Save problem features to find true synergy
	for w in workstations, t in times
		if (w in sp_workstations) & (t in sp_times)
			x_k[subproblem_id,w,t] = 1
			for p in impact_pods
				v_k[subproblem_id,p,w,t] = 1
			end
			for p in setdiff(pods, impact_pods)
				v_k[subproblem_id,p,w,t] = 0
			end
		else
			x_k[subproblem_id,w,t] = 0
			for p in pods
				v_k[subproblem_id,p,w,t] = 0
			end
		end
	end

	for m in orders, p in pods
		if (m in impact_orders) & (p in impact_pods)
			y_k[subproblem_id,m,p] = 1
		else
			y_k[subproblem_id,m,p] = 0
		end
	end

	for w in workstations, p in pods
		if (w in sp_workstations) & (p in impact_pods)
			z_k[subproblem_id,p,w] = 1
		else
			z_k[subproblem_id,p,w] = 0
		end
	end

	for m1 in orders, m2 in orders
		if (m1 in impact_orders) & (m2 in impact_orders)
			q_k[subproblem_id,m1,m2] = 1
		else
			q_k[subproblem_id,m1,m2] = 0
		end
	end

	obj_k[subproblem_id] = new_obj 

end