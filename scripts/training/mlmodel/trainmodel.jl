
include("models.jl") 

#-----------------------------------------------------------------------------------#

function initializedatastructures()

	w_features, t_features, m_features, p_features = Dict(), Dict(), Dict(), Dict()
	wt_features, mp_features, pw_features, pwt_features = Dict(), Dict(), Dict(), Dict(), Dict()
	actualobj, actualimp = Dict(), Dict()
	oldobj = Dict()
	sp_workstations, sp_pods, sp_times, sp_orders = Dict(), Dict(), Dict(), Dict(), Dict()

	return w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, actualobj, oldobj, actualimp, sp_workstations, sp_pods, sp_times, sp_orders

end

#-----------------------------------------------------------------------------------#

function featurebank()

	w_feat = ["w_centrality"] 
	t_feat = ["t_horizonrelativetime", "t_timeofday"]
	p_feat = ["p_avgdisttostations", "p_pctusefulitems"]
	m_feat = ["m_ordersize", "m_oneitemflag"]
	i_feat = ["i_iteminventory", "i_itemnumpods"]

	wt_feat = ["wt_queuecong", "wt_avglocalcong", "wt_maxlocalcong", "wt_avghyperlocalcong", "wt_maxhyperlocalcong"]
	mp_feat = ["mp_itemoverlap", "mp_itemoverlappct", "mp_overlapitemavginv", "mp_overlapitemavgaltpods"]
	pw_feat = ["pw_distance", "pw_xsteps", "pw_ysteps"]
	
	pwt_feat = ["pwt_existingcong"]

	featurenames = []
	for featgrp in [w_feat, t_feat, p_feat, m_feat, i_feat, wt_feat, mp_feat, pw_feat, pwt_feat]
		featurenames = union(featurenames, featgrp)
	end
	allfeatures = [j for j in 1:length(featurenames)]
	featurelookup = Dict()
	for f in 1:length(allfeatures)
		featurelookup[featurenames[f]] = f
	end

	w_featnums = [featurelookup[f] for f in w_feat]
	t_featnums = [featurelookup[f] for f in t_feat]
	p_featnums = [featurelookup[f] for f in p_feat]
	m_featnums = [featurelookup[f] for f in m_feat]
	wt_featnums = [featurelookup[f] for f in wt_feat]
	mp_featnums = [featurelookup[f] for f in mp_feat]
 	pw_featnums = [featurelookup[f] for f in pw_feat]
	pwt_featnums = [featurelookup[f] for f in pwt_feat]

	featuresforprediction_wt = union(wt_featnums, w_featnums, t_featnums)
	featuresforprediction_mp = union(mp_featnums, m_featnums, p_featnums)
	featuresforprediction_pw = union(pw_featnums, p_featnums, w_featnums)
	featuresforprediction_pwt = union(pwt_featnums, p_featnums, w_featnums, t_featnums)

	return allfeatures, featurenames, featurelookup, featuresforprediction_wt, featuresforprediction_mp, featuresforprediction_pw, featuresforprediction_pwt, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums

end

#-----------------------------------------------------------------------------------#

function addinstance(i, data, featurenames, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums)

	orders, workstations, pods, itemson, podswith, horizon, tstep = data["instance_features"]["orders"], data["instance_features"]["workstations"], data["instance_features"]["pods"], data["instance_features"]["itemson"], data["instance_features"]["podswith"], data["instance_features"]["horizon"], data["instance_features"]["tstep"]

	times = 0-tstep:tstep:horizon+tstep
	truetimes = 0:tstep:horizon

	for f in w_featnums
		for w in workstations
			w_features[i,f,w] = 0
		end
		for item in data["w_features"][replace(featurenames[f], "w_" => "" )]
			w = item[1]
			w_features[i,f,w] = item[2]
		end
	end

	for f in t_featnums
		for t in times
			t_features[i,f,t] = 0
		end
		for item in data["t_features"][replace(featurenames[f], "t_" => "" )]
			t = item[1]
			t_features[i,f,t] = item[2]
		end
	end

	for f in m_featnums
		for m in workstations
			m_features[i,f,m] = 0
		end
		for item in data["m_features"][replace(featurenames[f], "m_" => "" )]
			m = item[1]
			m_features[i,f,m] = item[2]
		end
	end

	for f in p_featnums
		for p in pods
			p_features[i,f,p] = 0
		end
		for item in data["p_features"][replace(featurenames[f], "p_" => "" )]
			p = item[1]
			p_features[i,f,p] = item[2]
		end
	end

	#-------------------------------------------#

	for f in wt_featnums
		for w in workstations, t in times
			wt_features[i,f,w,t] = 0
		end
		#for item in data["wt_features"][replace(featurenames[f], "wt_" => "" )]
		#	w, t = item[1]
		#	wt_features[i,f,w,t] = item[2]
		#end
		featuredata = data["wt_features"][replace(featurenames[f], "wt_" => "" )]
		for w in workstations, t in truetimes
			w_prime, t_prime = w - minimum(workstations) + 1, convert(Int, round(t/tstep + 1, digits=0))
			wt_features[i,f,w,t] = featuredata[w_prime, t_prime]
		end
	end

	for f in mp_featnums
		#for item in data["mp_features"][replace(featurenames[f], "mp_" => "" )]
		#	m, p = item[1]
		#	mp_features[i,f,m,p] = item[2]
		#end
		featuredata = data["mp_features"][replace(featurenames[f], "mp_" => "" )]
		for m in orders, p in pods
			mp_features[i,f,m,p] = featuredata[m, p]
		end
	end

	for f in pw_featnums
		#for item in data["pw_features"][replace(featurenames[f], "pw_" => "" )]
		#	p, w = item[1]
		#	pw_features[i,f,p,w] = item[2]
		#end
		featuredata = data["pw_features"][replace(featurenames[f], "pw_" => "" )]
		for p in pods, w in workstations
			w_prime = w - minimum(workstations) + 1
			pw_features[i,f,p,w] = featuredata[p, w_prime]
		end
	end

	for f in pwt_featnums
		for p in pods, w in workstations, t in times
			pwt_features[i,f,p,w,t] = 0
		end
		#for item in data["pwt_features"][replace(featurenames[f], "pwt_" => "" )]
		#	p, w, t = item[1]
		#	pwt_features[i,f,p,w,t] = item[2]
		#end
		featuredata = data["pwt_features"][replace(featurenames[f], "pwt_" => "" )]
		for p in pods, w in workstations, t in truetimes
			w_prime, t_prime = w - minimum(workstations) + 1, convert(Int, round(t/tstep + 1, digits=0))
			pwt_features[i,f,p,w,t] = featuredata[p, w_prime, t_prime]
		end
	end
	
	return w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features

end

#-----------------------------------------------------------------------------------#

function addsubproblem(i, k, spdata, tstep, targetnumpods, targetnumorders, sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp)
	
	sp_workstations[i,k] = []
	sp_times[i,k] = []
	sp_orders[i,k] = []
	sp_pods[i,k] = []

	for w in spdata["sp_workstations"]
		push!(sp_workstations[i,k], w)
	end

	impactorders, otherorders = spdata["impact_orders"], setdiff(spdata["sp_orders"], spdata["impact_orders"])
	numchosenorders = max(targetnumorders, length(impactorders))
	for m in spdata["impact_orders"]
		push!(sp_orders[i,k], m)
	end
	for mindex in randperm(length(otherorders))[1:min(length(otherorders), numchosenorders-length(impactorders))]
		push!(sp_orders[i,k], otherorders[mindex])
	end

	impactpods, otherpods = spdata["impact_pods"], setdiff(spdata["sp_pods"], spdata["impact_pods"])
	numchosenpods = max(targetnumpods, length(impactpods))
	for p in spdata["impact_pods"]
		push!(sp_pods[i,k], p)
	end
	for pindex in randperm(length(otherpods))[1:min(length(otherpods), numchosenpods-length(impactpods))]
		push!(sp_pods[i,k], otherpods[pindex])
	end

	for t in spdata["sp_tstart"]:tstep:spdata["sp_tend"]
		push!(sp_times[i,k], t)
	end

	actualobj[i,k] = spdata["new_obj"]
	oldobj[i,k] = spdata["old_obj"]
	actualimp[i,k] = spdata["new_obj"] - spdata["old_obj"]

	#mlpass[i,k] = spdata[]

	return sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp

end

#-----------------------------------------------------------------------------------#

function addinstance_nocongestion(i, data, featurenames, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums, featurelookup)

	println("Instance $i has no congestion")
	orders, workstations, pods, itemson, podswith, horizon, tstep = data["instance_features"]["orders"], data["instance_features"]["workstations"], data["instance_features"]["pods"], data["instance_features"]["itemson"], data["instance_features"]["podswith"], data["instance_features"]["horizon"], data["instance_features"]["tstep"]

	times = 0-tstep:tstep:horizon+tstep
	truetimes = 0:tstep:horizon

	for f in w_featnums
		for w in workstations
			w_features[i,f,w] = 0
		end
		for item in data["w_features"][replace(featurenames[f], "w_" => "" )]
			w = item[1]
			w_features[i,f,w] = item[2]
		end
	end

	for f in t_featnums
		for t in times
			t_features[i,f,t] = 0
		end
		for item in data["t_features"][replace(featurenames[f], "t_" => "" )]
			t = item[1]
			t_features[i,f,t] = item[2]
		end
	end

	for f in m_featnums
		for m in workstations
			m_features[i,f,m] = 0
		end
		for item in data["m_features"][replace(featurenames[f], "m_" => "" )]
			m = item[1]
			m_features[i,f,m] = item[2]
		end
	end

	for f in p_featnums
		for p in pods
			p_features[i,f,p] = 0
		end
		for item in data["p_features"][replace(featurenames[f], "p_" => "" )]
			p = item[1]
			p_features[i,f,p] = item[2]
		end
	end

	#-------------------------------------------#

	for f in wt_featnums
		for w in workstations, t in times
			wt_features[i,f,w,t] = 0
		end
		if !(f in [featurelookup["wt_queuecong"], featurelookup["wt_avglocalcong"], featurelookup["wt_maxlocalcong"], featurelookup["wt_maxhyperlocalcong"], featurelookup["wt_avghyperlocalcong"]])
			#for item in data["wt_features"][replace(featurenames[f], "wt_" => "" )]
			#	w, t = item[1]
			#	wt_features[i,f,w,t] = item[2]
			#end
			featuredata = data["wt_features"][replace(featurenames[f], "wt_" => "" )]
			for w in workstations, t in truetimes
				w_prime, t_prime = w - minimum(workstations) + 1, convert(Int, round(t/tstep + 1, digits=0))
				wt_features[i,f,w,t] = featuredata[w_prime, t_prime]
			end
		end
	end

	for f in mp_featnums
		#for item in data["mp_features"][replace(featurenames[f], "mp_" => "" )]
		#	m, p = item[1]
		#	mp_features[i,f,m,p] = item[2]
		#end
		featuredata = data["mp_features"][replace(featurenames[f], "mp_" => "" )]
		for m in orders, p in pods
			mp_features[i,f,m,p] = featuredata[m, p]
		end
	end

	for f in pw_featnums
		#for item in data["pw_features"][replace(featurenames[f], "pw_" => "" )]
		#	p, w = item[1]
		#	pw_features[i,f,p,w] = item[2]
		#end
		featuredata = data["pw_features"][replace(featurenames[f], "pw_" => "" )]
		for p in pods, w in workstations
			w_prime = w - minimum(workstations) + 1
			pw_features[i,f,p,w] = featuredata[p, w_prime]
		end
	end

	for f in pwt_featnums
		for p in pods, w in workstations, t in times
			pwt_features[i,f,p,w,t] = 0
		end
		if f != featurelookup["pwt_existingcong"]
			#for item in data["pwt_features"][replace(featurenames[f], "pwt_" => "" )]
			#	p, w, t = item[1]
			#	pwt_features[i,f,p,w,t] = item[2]
			#end
			featuredata = data["pwt_features"][replace(featurenames[f], "pwt_" => "" )]
			for p in pods, w in workstations, t in truetimes
				w_prime, t_prime = w - minimum(workstations) + 1, convert(Int, round(t/tstep + 1, digits=0))
				pwt_features[i,f,p,w,t] = featuredata[p, w_prime, t_prime]
			end
		end
	end
	
	return w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features

end

#-----------------------------------------------------------------------------------#

function addsubproblem_nocongestion(i, k, spdata, tstep, targetnumpods, targetnumorders, sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp)
	
	sp_workstations[i,k] = []
	sp_times[i,k] = []
	sp_orders[i,k] = []
	sp_pods[i,k] = []

	for w in spdata["sp_workstations"]
		push!(sp_workstations[i,k], w)
	end

	impactorders, otherorders = spdata["nocongestion_impact_orders"], setdiff(spdata["sp_orders"], spdata["nocongestion_impact_orders"])
	numchosenorders = max(targetnumorders, length(impactorders))
	for m in spdata["nocongestion_impact_orders"]
		push!(sp_orders[i,k], m)
	end
	for mindex in randperm(length(otherorders))[1:min(length(otherorders), numchosenorders-length(impactorders))]
		push!(sp_orders[i,k], otherorders[mindex])
	end

	impactpods, otherpods = spdata["nocongestion_impact_pods"], setdiff(spdata["sp_pods"], spdata["nocongestion_impact_pods"])
	numchosenpods = max(targetnumpods, length(impactpods))
	for p in spdata["nocongestion_impact_pods"]
		push!(sp_pods[i,k], p)
	end
	for pindex in randperm(length(otherpods))[1:min(length(otherpods), numchosenpods-length(impactpods))]
		push!(sp_pods[i,k], otherpods[pindex])
	end

	for t in spdata["sp_tstart"]:tstep:spdata["sp_tend"]
		push!(sp_times[i,k], t)
	end

	actualobj[i,k] = spdata["nocongestion_obj"]
	oldobj[i,k] = spdata["old_obj"]
	actualimp[i,k] = spdata["nocongestion_obj"] - spdata["old_obj"]

	return sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp

end

#-----------------------------------------------------------------------------------#

function findbetaestimates(sp_workstations, sp_pods, sp_orders, sp_times, actualobj)

	#consistency_flag, x, y, z, v, q, actualobj = 1, length(workstationgrouplist), x_k, y_k, z_k, v_k, q_k, obj_k

	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", 60*60) # iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", 1)

	#Variables
	@variable(model, pred[i = probleminstances, k = 1:numsubproblems[i]])
	@variable(model, err[i = probleminstances, k = 1:numsubproblems[i]])
	@variable(model, beta_wt[f = featuresforprediction_wt])
	@variable(model, beta_mp[f = featuresforprediction_mp])
	@variable(model, beta_pw[f = featuresforprediction_pw])
	@variable(model, beta_pwt[f = featuresforprediction_pwt])

	@variable(model, nonzero_wt[f = featuresforprediction_wt], Bin)
	@variable(model, nonzero_mp[f = featuresforprediction_mp], Bin)
	@variable(model, nonzero_pw[f = featuresforprediction_pw], Bin)
	@variable(model, nonzero_pwt[f = featuresforprediction_pwt], Bin)
	@variable(model, chooseone, Bin)

	#Objective
	@objective(model, Min, sum(sum(err[i,k] for k in 1:numsubproblems[i]) for i in probleminstances))

	#Prediction calculation
	@constraint(model, predictioncalc[i = probleminstances, k = 1:numsubproblems[i]], pred[i,k] == 
		#Workstation-time
		sum(sum(sum(beta_wt[f] * wt_features[i,f,w,t] for f in wt_featnums) for w in sp_workstations[i,k]) for t in sp_times[i,k]) 
		+ sum(sum(beta_wt[f] * w_features[i,f,w] for f in w_featnums) for w in sp_workstations[i,k]) 
		+ sum(sum(beta_wt[f] * t_features[i,f,t] for f in t_featnums) for t in sp_times[i,k]) 
		#Order-pod
		+ sum(sum(sum(beta_mp[f] * mp_features[i,f,m,p] for f in mp_featnums) for m in sp_orders[i,k]) for p in sp_pods[i,k]) 
		+ sum(sum(beta_mp[f] * m_features[i,f,m] for f in m_featnums) for m in sp_orders[i,k]) 
		+ sum(sum(beta_mp[f] * p_features[i,f,p] for f in p_featnums) for p in sp_pods[i,k]) 
		#Pod-workstation
		+ sum(sum(sum(beta_pw[f] * pw_features[i,f,p,w] for f in pw_featnums) for w in sp_workstations[i,k]) for p in sp_pods[i,k]) 
		+ sum(sum(beta_pw[f] * p_features[i,f,p] for f in p_featnums) for p in sp_pods[i,k]) 
		+ sum(sum(beta_pw[f] * w_features[i,f,w] for f in w_featnums) for w in sp_workstations[i,k]) 
		#Pod-workstation-time
		+ sum(sum(sum(sum(beta_pwt[f] * pwt_features[i,f,p,w,t] for f in pwt_featnums) for w in sp_workstations[i,k]) for p in sp_pods[i,k]) for t in sp_times[i,k])  
		+ sum(sum(beta_pwt[f] * p_features[i,f,p] for f in p_featnums) for p in sp_pods[i,k])
		+ sum(sum(beta_pwt[f] * w_features[i,f,w] for f in w_featnums) for w in sp_workstations[i,k])   
		+ sum(sum(beta_pwt[f] * t_features[i,f,t] for f in t_featnums) for t in sp_times[i,k])  
		)

	#for i in probleminstances, k in 1:numsubproblems[i]
	#	println(i, ", ", k)
	#	@constraint(model, pred[i,k] == sum(sum(sum(beta_wt[f] * wt_features[i,f,w,t] for f in wt_featnums) for w in sp_workstations[i,k]) for t in sp_times[i,k]) + sum(sum(beta_wt[f] * w_features[i,f,w] for f in w_featnums) for w in sp_workstations[i,k]) + sum(sum(beta_wt[f] * t_features[i,f,t] for f in t_featnums) for t in sp_times[i,k]) + sum(sum(sum(beta_mp[f] * mp_features[i,f,m,p] for f in mp_featnums) for m in sp_orders[i,k]) for p in sp_pods[i,k]) + sum(sum(beta_mp[f] * m_features[i,f,m] for f in m_featnums) for m in sp_orders[i,k]) + sum(sum(beta_mp[f] * p_features[i,f,p] for f in p_featnums) for p in sp_pods[i,k]) + sum(sum(sum(beta_pw[f] * pw_features[i,f,p,w] for f in pw_featnums) for w in sp_workstations[i,k]) for p in sp_pods[i,k]) + sum(sum(beta_pw[f] * p_features[i,f,p] for f in p_featnums) for p in sp_pods[i,k]) + sum(sum(beta_pw[f] * w_features[i,f,w] for f in w_featnums) for w in sp_workstations[i,k]) + sum(sum(sum(sum(beta_pwt[f] * pwt_features[i,f,p,w,t] for f in pwt_featnums) for w in sp_workstations[i,k]) for p in sp_pods[i,k]) for t in sp_times[i,k])  + sum(sum(beta_pwt[f] * p_features[i,f,p] for f in p_featnums) for p in sp_pods[i,k]) + sum(sum(beta_pwt[f] * w_features[i,f,w] for f in w_featnums) for w in sp_workstations[i,k]) + sum(sum(beta_pwt[f] * t_features[i,f,t] for f in t_featnums) for t in sp_times[i,k])  )
	#end

	#Prediction error
	@constraint(model, prederror1[i = probleminstances, k = 1:numsubproblems[i]], err[i,k] >= pred[i,k] - actualobj[i,k])
	@constraint(model, prederror2[i = probleminstances, k = 1:numsubproblems[i]], err[i,k] >= actualobj[i,k] - pred[i,k])

	#Regularization

	#Consistency / interpretability of predictions
	#@constraint(model, otherinstancesareworse[i = probleminstances, k = 1:numsubproblems[i]], pred[i, k] <= pred[i, k_star])
	#@constraint(model, nonzerofeatures_wt_1[f in featuresforprediction_wt], beta_wt[f] <= 10000 * nonzero_wt[f])
	#@constraint(model, nonzerofeatures_wt_2[f in featuresforprediction_wt], beta_wt[f] >= -10000 * nonzero_wt[f])
	@constraint(model, nonzerofeatures_mp_1[f in featuresforprediction_mp], beta_mp[f] <= 10000 * nonzero_mp[f])
	@constraint(model, nonzerofeatures_mp_2[f in featuresforprediction_mp], beta_mp[f] >= -10000 * nonzero_mp[f])
	#@constraint(model, nonzerofeatures_pw_1[f in featuresforprediction_pw], beta_pw[f] <= 10000 * nonzero_pw[f])
	#@constraint(model, nonzerofeatures_pw_2[f in featuresforprediction_pw], beta_pw[f] >= -10000 * nonzero_pw[f])
	#@constraint(model, nonzerofeatures_pwt_1[f in featuresforprediction_pwt], beta_pwt[f] <= 10000 * nonzero_pwt[f])
	#@constraint(model, nonzerofeatures_pwt_2[f in featuresforprediction_pwt], beta_pwt[f] >= -10000 * nonzero_pwt[f])
	@constraint(model, nonzero_mp[featurelookup["mp_itemoverlap"]] + nonzero_mp[featurelookup["mp_itemoverlappct"]] <= 1)
	@constraint(model, nonzero_mp[featurelookup["m_ordersize"]] + nonzero_mp[featurelookup["m_oneitemflag"]] <= 1)
	@constraint(model, beta_pw[featurelookup["pw_itemoverlap"]] == 0)

	#====================================================#

	#beta_wt, beta_mp, beta_pw, beta_pwt = value.(beta_wt), value.(beta_mp), value.(beta_pw), value.(beta_pwt)

	#Solve IP
	status = optimize!(model)

	if termination_status(model) == MOI.OPTIMAL
		solvetime = solve_time(model)
		obj = objective_value(model)
		feasibleflag = 1
		#println("Total items picked = ", obj)
	else
		println("Error in solving!")
		solvetime = solve_time(model)
		obj = 0
		#println("Total items picked = ", obj)
		feasibleflag = 0
	end

	#====================================================#

	return value.(beta_wt), value.(beta_mp), value.(beta_pw), value.(beta_pwt), value.(pred)

end

#-----------------------------------------------------------------------------------#

function calcpredictionfortestset(testinstances, beta_wt, beta_mp, beta_pw, beta_pwt, numsubproblems, sp_workstations, sp_times, sp_pods, sp_orders, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums)
 
	prediction = Dict() 

	for i in testinstances, k in 1:numsubproblems[i]

		prediction[i,k] = 0

		for t in sp_times[i,k], w in sp_workstations[i,k], f in wt_featnums
			prediction[i,k] += beta_wt[f] * wt_features[i,f,w,t]
		end
		for w in sp_workstations[i,k], f in w_featnums
			prediction[i,k] += beta_wt[f] * w_features[i,f,w]
		end
		for t in sp_times[i,k], f in t_featnums
			prediction[i,k] += beta_wt[f] * t_features[i,f,t]
		end

		for p in sp_pods[i,k], m in sp_orders[i,k], f in mp_featnums
			prediction[i,k] += beta_mp[f] * mp_features[i,f,m,p]
		end
		for m in sp_orders[i,k], f in m_featnums
			prediction[i,k] += beta_mp[f] * m_features[i,f,m]
		end
		for p in sp_pods[i,k], f in p_featnums
			prediction[i,k] += beta_mp[f] * p_features[i,f,p]
		end

		for p in sp_pods[i,k],  w in sp_workstations[i,k], f in pw_featnums
			prediction[i,k] += beta_pw[f] * pw_features[i,f,p,w]
		end
		for w in sp_workstations[i,k], f in w_featnums
			prediction[i,k] += beta_pw[f] * w_features[i,f,w]
		end
		for p in sp_pods[i,k], f in p_featnums
			prediction[i,k] += beta_pw[f] * p_features[i,f,p]
		end

		for p in sp_pods[i,k], w in sp_workstations[i,k], t in sp_times[i,k], f in pwt_featnums
			prediction[i,k] += beta_pwt[f] * pwt_features[i,f,p,w,t]
		end
		for p in sp_pods[i,k], f in p_featnums
			prediction[i,k] += beta_pwt[f] * p_features[i,f,p]
		end
		for w in sp_workstations[i,k], f in w_featnums
			prediction[i,k] += beta_pwt[f] * w_features[i,f,w]
		end
		for t in sp_times[i,k], f in t_featnums
			prediction[i,k] += beta_pwt[f] * t_features[i,f,t]
		end

	end

	return prediction

end

#-----------------------------------------------------------------------------------#

# testinstances, beta_wt, beta_mp, beta_pw, beta_pwt, shifts, numsubproblems, sp_workstations, sp_times, sp_pods, sp_orders, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums = testinstances, mlmodel.beta_wt, mlmodel.beta_mp, mlmodel.beta_pw, mlmodel.beta_pwt, mlmodel.shifts, numsubproblems, sp_workstations, sp_times, sp_pods, sp_orders, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums
#for f in 1:numfeatures
#	println(featurenames[f], " --> ", featurevalues[i,f], " * ", beta[f], " = ", featurevalues[i,f] * beta[f])
#end  
#=

if 1==1
for f in wt_featnums
	featval = 0
	totalval = 0
	for t in sp_times[i,k], w in sp_workstations[i,k]
		featval +=  wt_features[i,f,w,t]
		totalval += beta_wt[f] * wt_features[i,f,w,t]
	end
	println(featurenames[f], " --> ", featval, " * ",  beta_wt[f], " = ", totalval)
end
for f in w_featnums
	featval = 0
	totalval = 0
	for w in sp_workstations[i,k]
		featval +=  w_features[i,f,w]
		totalval += beta_wt[f] * w_features[i,f,w]
	end
	println(featurenames[f], " --> ", featval, " * ",  beta_wt[f], " = ", totalval)
end
for f in t_featnums
	featval = 0
	totalval = 0
	for t in sp_times[i,k]
		featval +=  t_features[i,f,t]
		totalval += beta_wt[f] * t_features[i,f,t]
	end
	println(featurenames[f], " --> ", featval, " * ",  beta_wt[f], " = ", totalval)
end
println("--------------------------")
for f in mp_featnums
	featval = 0
	totalval = 0
	for p in sp_pods[i,k], m in sp_orders[i,k]
		featval +=  mp_features[i,f,m,p]
		totalval += beta_mp[f] * mp_features[i,f,m,p]
	end
	println(featurenames[f], " --> ", featval, " * ",  beta_mp[f], " = ", totalval)
end
for f in m_featnums
	featval = 0
	totalval = 0
	for m in sp_orders[i,k]
		featval +=  m_features[i,f,m]
		totalval += beta_mp[f] * m_features[i,f,m]
	end
	println(featurenames[f], " --> ", featval, " * ",  beta_mp[f], " = ", totalval)
end
for f in p_featnums
	featval = 0
	totalval = 0
	for p in sp_pods[i,k]
		featval +=  p_features[i,f,p]
		totalval += beta_mp[f] * p_features[i,f,p]
	end
	println(featurenames[f], " --> ", featval, " * ",  beta_mp[f], " = ", totalval)
end
println("--------------------------")
for f in pw_featnums
	featval = 0
	totalval = 0
	for p in sp_pods[i,k], w in sp_workstations[i,k]
		featval +=  pw_features[i,f,p,w]
		totalval += beta_pw[f] * pw_features[i,f,p,w]
	end
	println(featurenames[f], " --> ", featval, " * ",  beta_pw[f], " = ", totalval)
end
for f in w_featnums
	featval = 0
	totalval = 0
	for w in sp_workstations[i,k]
		featval +=  w_features[i,f,w]
		totalval += beta_pw[f] * w_features[i,f,w]
	end
	println(featurenames[f], " --> ", featval, " * ",  beta_pw[f], " = ", totalval)
end
for f in p_featnums
	featval = 0
	totalval = 0
	for p in sp_pods[i,k]
		featval +=  p_features[i,f,p]
		totalval += beta_pw[f] * p_features[i,f,p]
	end
	println(featurenames[f], " --> ", featval, " * ",  beta_pw[f], " = ", totalval)
end
println("--------------------------")
for f in pwt_featnums
	featval = 0
	totalval = 0
	for p in sp_pods[i,k], w in sp_workstations[i,k], t in sp_times[i,k]
		featval +=  pwt_features[i,f,p,w, t]
		totalval += beta_pwt[f] * pwt_features[i,f,p,w,t]
	end
	println(featurenames[f], " --> ", featval, " * ",  beta_pwt[f], " = ", totalval)
end
for f in w_featnums
	featval = 0
	totalval = 0
	for w in sp_workstations[i,k]
		featval +=  w_features[i,f,w]
		totalval += beta_pwt[f] * w_features[i,f,w]
	end
	println(featurenames[f], " --> ", featval, " * ",  beta_pwt[f], " = ", totalval)
end
for f in p_featnums
	featval = 0
	totalval = 0
	for p in sp_pods[i,k]
		featval +=  p_features[i,f,p]
		totalval += beta_pwt[f] * p_features[i,f,p]
	end
	println(featurenames[f], " --> ", featval, " * ",  beta_pwt[f], " = ", totalval)
end
for f in t_featnums
	featval = 0
	totalval = 0
	for t in sp_times[i,k]
		featval +=  t_features[i,f,t]
		totalval += beta_pwt[f] * t_features[i,f,t]
	end
	println(featurenames[f], " --> ", featval, " * ",  beta_pwt[f], " = ", totalval)
end
println("--------------------------")
for w in sp_workstations[i,k]
	println("w_shift --> ", shifts[3])
end
for t in sp_times[i,k]
	println("t_shift --> ", shifts[4])
end
end
=#

function calcpredictionfortestset_shift(testinstances, beta_wt, beta_mp, beta_pw, beta_pwt, shifts, numsubproblems, sp_workstations, sp_times, sp_pods, sp_orders, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums)
 
	prediction = Dict() 

	for i in testinstances, k in 1:numsubproblems[i]

		prediction[i,k] = 0

		for t in sp_orders[i,k]
			prediction[i,k] += shifts[1]
		end
		for t in sp_pods[i,k]
			prediction[i,k] += shifts[2]
		end
		for t in sp_workstations[i,k]
			prediction[i,k] += shifts[3]
		end
		for t in sp_times[i,k]
			prediction[i,k] += shifts[4]
		end

		for t in sp_times[i,k], w in sp_workstations[i,k], f in wt_featnums
			prediction[i,k] += beta_wt[f] * wt_features[i,f,w,t]
		end
		for w in sp_workstations[i,k], f in w_featnums
			prediction[i,k] += beta_wt[f] * w_features[i,f,w]
		end
		for t in sp_times[i,k], f in t_featnums
			prediction[i,k] += beta_wt[f] * t_features[i,f,t]
		end

		for p in sp_pods[i,k], m in sp_orders[i,k], f in mp_featnums
			prediction[i,k] += beta_mp[f] * mp_features[i,f,m,p]
		end
		for m in sp_orders[i,k], f in m_featnums
			prediction[i,k] += beta_mp[f] * m_features[i,f,m]
		end
		for p in sp_pods[i,k], f in p_featnums
			prediction[i,k] += beta_mp[f] * p_features[i,f,p]
		end

		for p in sp_pods[i,k],  w in sp_workstations[i,k], f in pw_featnums
			prediction[i,k] += beta_pw[f] * pw_features[i,f,p,w]
		end
		for w in sp_workstations[i,k], f in w_featnums
			prediction[i,k] += beta_pw[f] * w_features[i,f,w]
		end
		for p in sp_pods[i,k], f in p_featnums
			prediction[i,k] += beta_pw[f] * p_features[i,f,p]
		end

		for p in sp_pods[i,k], w in sp_workstations[i,k], t in sp_times[i,k], f in pwt_featnums
			prediction[i,k] += beta_pwt[f] * pwt_features[i,f,p,w,t]
		end
		for p in sp_pods[i,k], f in p_featnums
			prediction[i,k] += beta_pwt[f] * p_features[i,f,p]
		end
		for w in sp_workstations[i,k], f in w_featnums
			prediction[i,k] += beta_pwt[f] * w_features[i,f,w]
		end
		for t in sp_times[i,k], f in t_featnums
			prediction[i,k] += beta_pwt[f] * t_features[i,f,t]
		end

	end

	return prediction

end

#-----------------------------------------------------------------------------------#

function calcpredictionfortestset_compare(testinstances, beta_wt, beta_mp, beta_pw, beta_pwt, shifts, alpha_m, alpha_p, alpha_w, alpha_t, numsubproblems, sp_workstations, sp_times, sp_pods, sp_orders, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums)
 
	prediction = Dict() 

	for i in testinstances, k in 1:numsubproblems[i]

		prediction[i,k] = 0

		for t in sp_orders[i,k]
			prediction[i,k] += shifts[1]
		end
		for t in sp_pods[i,k]
			prediction[i,k] += shifts[2]
		end
		for t in sp_workstations[i,k]
			prediction[i,k] += shifts[3]
		end
		for t in sp_times[i,k]
			prediction[i,k] += shifts[4]
		end

		for t in sp_times[i,k], w in sp_workstations[i,k], f in wt_featnums
			prediction[i,k] += beta_wt[f] * wt_features[i,f,w,t]
		end
		for w in sp_workstations[i,k], f in w_featnums
			prediction[i,k] += beta_wt[f] * w_features[i,f,w]
		end
		for t in sp_times[i,k], f in t_featnums
			prediction[i,k] += beta_wt[f] * t_features[i,f,t]
		end

		for p in sp_pods[i,k], m in sp_orders[i,k], f in mp_featnums
			prediction[i,k] += beta_mp[f] * mp_features[i,f,m,p]
		end
		for m in sp_orders[i,k], f in m_featnums
			prediction[i,k] += beta_mp[f] * m_features[i,f,m]
		end
		for p in sp_pods[i,k], f in p_featnums
			prediction[i,k] += beta_mp[f] * p_features[i,f,p]
		end

		for p in sp_pods[i,k],  w in sp_workstations[i,k], f in pw_featnums
			prediction[i,k] += beta_pw[f] * pw_features[i,f,p,w]
		end
		for w in sp_workstations[i,k], f in w_featnums
			prediction[i,k] += beta_pw[f] * w_features[i,f,w]
		end
		for p in sp_pods[i,k], f in p_featnums
			prediction[i,k] += beta_pw[f] * p_features[i,f,p]
		end

		for p in sp_pods[i,k], w in sp_workstations[i,k], t in sp_times[i,k], f in pwt_featnums
			prediction[i,k] += beta_pwt[f] * pwt_features[i,f,p,w,t]
		end
		for p in sp_pods[i,k], f in p_featnums
			prediction[i,k] += beta_pwt[f] * p_features[i,f,p]
		end
		for w in sp_workstations[i,k], f in w_featnums
			prediction[i,k] += beta_pwt[f] * w_features[i,f,w]
		end
		for t in sp_times[i,k], f in t_featnums
			prediction[i,k] += beta_pwt[f] * t_features[i,f,t]
		end

	end

	return prediction

end

#-----------------------------------------------------------------------------------#

function traindynamicsubproblemmodel(trainfilelist, traincongestion_flag, trainassignment_flag, setbetas)

	w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, actualobj, oldobj, actualimp, sp_workstations, sp_times, sp_pods, sp_orders = initializedatastructures()

	probleminstances = []
	numsubproblems, weight = Dict(), Dict()

	#Add data for each instance
	for file in trainfilelist #Each file contains one instance and some number of subproblems 

		println(file)

		datafile = string(trainingfolder, file)
		data = load(datafile)

		i = parse(Int64, file[last(findfirst("_run", file))+1 : findfirst(".jld2", file)[1]-1 ])

		push!(probleminstances, i)
		push!(probleminstances, i+1000)
		numsubproblems[i] = length(data["subproblemdata"])
		numsubproblems[i+1000] = length(data["subproblemdata"])
		mlpass = 0

		#Add instance features
		w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features = addinstance(i, data, featurenames, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums)
		w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features = addinstance_nocongestion(i+1000, data, featurenames, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums, featurelookup)
		
		#Add subproblem structure and objectives
		for k in 1:numsubproblems[i]
			spdata = data["subproblemdata"][k]
			sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp = addsubproblem(i, k, spdata, tstep, targetnumpods, targetnumorders, sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp)
			weight[i,k] = 1
			sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp = addsubproblem_nocongestion(i+1000, k, spdata, tstep, targetnumpods, targetnumorders, sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp)
			weight[i+1000,k] = 1
		end

	end

	#---------------------------------------------------------------------------------------------------------#

	include("scripts/training/mlmodel/models.jl")

	#---------------------------------------------------------------------------------------------------------#

	#Train model without any congestion parameters
	beta_wt, beta_mp, beta_pw, beta_pwt, pred, shifts = linearregressionmodel_nocongestion(sp_workstations, sp_pods, sp_orders, sp_times, actualobj, probleminstances, numsubproblems, featurelookup, featuresforprediction_wt, featuresforprediction_mp, featuresforprediction_pw, featuresforprediction_pwt, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums, setbetas, traincongestion_flag, trainassignment_flag, weight)
	#printmodelsummary(beta_wt, beta_mp, beta_pw, beta_pwt, pred)

	println("------------------------------------------")
	for f in featuresforprediction_wt
		println(featurenames[f], " = ", value(beta_wt[f]))
	end
	println("------------------------------------------")
	for f in featuresforprediction_mp
		println(featurenames[f], " = ", value(beta_mp[f]))
	end
	println("------------------------------------------")
	for f in featuresforprediction_pw
		println(featurenames[f], " = ", value(beta_pw[f]))
	end
	println("------------------------------------------")
	for f in featuresforprediction_pwt
		println(featurenames[f], " = ", value(beta_pwt[f]))
	end
	println("------------------------------------------")
	println("beta_m = ", shifts[1])
	println("beta_p = ", shifts[2])
	println("beta_w = ", shifts[3])
	println("beta_t = ", shifts[4])
	println("------------------------------------------")

	totalAE, totalSE, tss, validsps = 0, 0, 0, 0
	meanobj = mean(values(actualobj))
	for i in probleminstances, k in 1:numsubproblems[i]
		if actualobj[i,k] > 1e-4
			totalAE += abs(pred[i,k] - actualobj[i,k])
			totalSE += (pred[i,k] - actualobj[i,k])^2
			tss += (meanobj - actualobj[i,k])^2
			validsps += 1
		end
	end

	println("Training mean absolute error = ", totalAE / validsps)
	println("Training mean squared error = ", totalSE / validsps)
	println("Training R_squared = ", 1 - totalSE/tss)

	mlmodel = (features = allfeatures, featurenames = featurenames, beta_wt = beta_wt, beta_mp = beta_mp, beta_pw = beta_pw, beta_pwt = beta_pwt, shifts = shifts, trainingpredictions = pred) 

	#---------------------------------------------------------------------------------------------------------#

	#Train model with congestion parameters
	beta_wt2, beta_mp2, beta_pw2, beta_pwt2, pred2, shifts2 = linearregressionmodel_wt(sp_workstations, sp_pods, sp_orders, sp_times, actualobj, probleminstances, numsubproblems, featurelookup, featuresforprediction_wt, featuresforprediction_mp, featuresforprediction_pw, featuresforprediction_pwt, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums, setbetas, traincongestion_flag, trainassignment_flag, weight)
	#printmodelsummary(beta_wt, beta_mp, beta_pw, beta_pwt, pred)

	println("------------------------------------------")
	for f in featuresforprediction_wt
		println(featurenames[f], " = ", value(beta_wt2[f]))
	end
	println("------------------------------------------")
	for f in featuresforprediction_mp
		println(featurenames[f], " = ", value(beta_mp2[f]))
	end
	println("------------------------------------------")
	for f in featuresforprediction_pw
		println(featurenames[f], " = ", value(beta_pw2[f]))
	end
	println("------------------------------------------")
	for f in featuresforprediction_pwt
		println(featurenames[f], " = ", value(beta_pwt2[f]))
	end
	println("------------------------------------------")
	println("beta_m = ", shifts2[1])
	println("beta_p = ", shifts2[2])
	println("beta_w = ", shifts2[3])
	println("beta_t = ", shifts2[4])
	println("------------------------------------------")

	totalAE, totalSE, tss, validsps = 0, 0, 0, 0
	meanobj = mean(values(actualobj))
	for i in probleminstances, k in 1:numsubproblems[i]
		if actualobj[i,k] > 1e-4
			totalAE += abs(pred2[i,k] - actualobj[i,k])
			totalSE += (pred2[i,k] - actualobj[i,k])^2
			tss += (meanobj - actualobj[i,k])^2
			validsps += 1
		end
	end

	println("Training mean absolute error = ", totalAE / validsps)
	println("Training mean squared error = ", totalSE / validsps)
	println("Training R_squared = ", 1 - totalSE/tss)

	mlmodel2 = (features = allfeatures, featurenames = featurenames, beta_wt = beta_wt2, beta_mp = beta_mp2, beta_pw = beta_pw2, beta_pwt = beta_pwt2, shifts = shifts2, trainingpredictions = pred2) 

	#---------------------------------------------------------------------------------------------------------#

	#Train model with intercept term (and congestion parameters)
	beta_wt3, beta_mp3, beta_pw3, beta_pwt3, pred3, shifts3 = linearregressionmodel_intercept(sp_workstations, sp_pods, sp_orders, sp_times, actualobj, probleminstances, numsubproblems, featurelookup, featuresforprediction_wt, featuresforprediction_mp, featuresforprediction_pw, featuresforprediction_pwt, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums, setbetas, traincongestion_flag, trainassignment_flag, weight)
	#printmodelsummary(beta_wt, beta_mp, beta_pw, beta_pwt, pred)

	println("------------------------------------------")
	for f in featuresforprediction_wt
		println(featurenames[f], " = ", value(beta_wt3[f]))
	end
	println("------------------------------------------")
	for f in featuresforprediction_mp
		println(featurenames[f], " = ", value(beta_mp3[f]))
	end
	println("------------------------------------------")
	for f in featuresforprediction_pw
		println(featurenames[f], " = ", value(beta_pw3[f]))
	end
	println("------------------------------------------")
	for f in featuresforprediction_pwt
		println(featurenames[f], " = ", value(beta_pwt3[f]))
	end
	println("------------------------------------------")
	println("beta_m = ", shifts3[1])
	println("beta_p = ", shifts3[2])
	println("beta_w = ", shifts3[3])
	println("beta_t = ", shifts3[4])
	println("beta_0 = ", shifts3[5])
	println("------------------------------------------")

	totalAE, totalSE, tss, validsps = 0, 0, 0, 0
	meanobj = mean(values(actualobj))
	for i in probleminstances, k in 1:numsubproblems[i]
		if actualobj[i,k] > 1e-4
			totalAE += abs(pred3[i,k] - actualobj[i,k])
			totalSE += (pred3[i,k] - actualobj[i,k])^2
			tss += (meanobj - actualobj[i,k])^2
			validsps += 1
		end
	end

	println("Training mean absolute error = ", totalAE / validsps)
	println("Training mean squared error = ", totalSE / validsps)
	println("Training R_squared = ", 1 - totalSE/tss)

	mlmodel3 = (features = allfeatures, featurenames = featurenames, beta_wt = beta_wt3, beta_mp = beta_mp3, beta_pw = beta_pw3, beta_pwt = beta_pwt3, shifts = shifts3, trainingpredictions = pred3) 

	#---------------------------------------------------------------------------------------------------------#

	return mlmodel, mlmodel2, mlmodel3

end

#-----------------------------------------------------------------------------------#

function traindynamicsubproblemmodel_compare(trainfilelist, traincongestion_flag, trainassignment_flag, setbetas)

	w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, actualobj, oldobj, actualimp, sp_workstations, sp_times, sp_pods, sp_orders = initializedatastructures()

	probleminstances = []
	numsubproblems, weight = Dict(), Dict()

	#Add data for each instance
	for file in trainfilelist #Each file contains one instance and some number of subproblems 

		println(file)

		datafile = string(trainingfolder, file)
		data = load(datafile)

		i = parse(Int64, file[last(findfirst("_run", file))+1 : findfirst(".jld2", file)[1]-1 ])

		push!(probleminstances, i)
		#push!(probleminstances, i+1000)
		numsubproblems[i] = length(data["subproblemdata"])
		numsubproblems[i+1000] = length(data["subproblemdata"])
		mlpass = 0

		#Add instance features
		w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features = addinstance(i, data, featurenames, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums)
		w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features = addinstance_nocongestion(i+1000, data, featurenames, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums, featurelookup)
		
		#Add subproblem structure and objectives
		for k in 1:numsubproblems[i]
			spdata = data["subproblemdata"][k]
			sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp = addsubproblem(i, k, spdata, tstep, targetnumpods, targetnumorders, sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp)
			weight[i,k] = 1
			sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp = addsubproblem_nocongestion(i+1000, k, spdata, tstep, targetnumpods, targetnumorders, sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp)
			weight[i+1000,k] = 1
		end

	end

	include("scripts/dynamicsubproblem/mlmodels.jl")
	beta_wt, beta_mp, beta_pw, beta_pwt, pred, shifts, alpha_m, alpha_p, alpha_w, alpha_t = linearregressionmodel_compare(sp_workstations, sp_pods, sp_orders, sp_times, actualobj, probleminstances, numsubproblems, featurelookup, featuresforprediction_wt, featuresforprediction_mp, featuresforprediction_pw, featuresforprediction_pwt, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums, setbetas, traincongestion_flag, trainassignment_flag, weight)
	#printmodelsummary(beta_wt, beta_mp, beta_pw, beta_pwt, pred)

	println("------------------------------------------")
	for f in wt_featnums
		println(featurenames[f], " = ", value(beta_wt[f]))
	end
	println("------------------------------------------")
	for f in mp_featnums
		println(featurenames[f], " = ", value(beta_mp[f]))
	end
	println("------------------------------------------")
	for f in pw_featnums
		println(featurenames[f], " = ", value(beta_pw[f]))
	end
	println("------------------------------------------")
	for f in pwt_featnums
		println(featurenames[f], " = ", value(beta_pwt[f]))
	end
	println("------------------------------------------")
	for f in m_featnums
		println(featurenames[f], " = ", value(alpha_m[f]))
	end
	println("------------------------------------------")
	for f in p_featnums
		println(featurenames[f], " = ", value(alpha_p[f]))
	end
	println("------------------------------------------")
	for f in w_featnums
		println(featurenames[f], " = ", value(alpha_w[f]))
	end
	println("------------------------------------------")
	for f in t_featnums
		println(featurenames[f], " = ", value(alpha_t[f]))
	end
	println("------------------------------------------")
	println("beta_m = ", shifts[1])
	println("beta_p = ", shifts[2])
	println("beta_w = ", shifts[3])
	println("beta_t = ", shifts[4])
	println("------------------------------------------")

	totalAE, totalSE, tss, validsps = 0, 0, 0, 0
	meanobj = mean(values(actualobj))
	for i in probleminstances, k in 1:numsubproblems[i]
		if actualobj[i,k] > 1e-4
			totalAE += abs(pred[i,k] - actualobj[i,k])
			totalSE += (pred[i,k] - actualobj[i,k])^2
			tss += (meanobj - actualobj[i,k])^2
			validsps += 1
		end
	end

	println("Training mean absolute error = ", totalAE / validsps)
	println("Training mean squared error = ", totalSE / validsps)
	println("Training R_squared = ", 1 - totalSE/tss)

	mlmodel = (features = allfeatures, featurenames = featurenames, beta_wt = beta_wt, beta_mp = beta_mp, beta_pw = beta_pw, beta_pwt = beta_pwt, shifts = shifts, alpha_m = alpha_m, alpha_p = alpha_p, alpha_w = alpha_w, alpha_t = alpha_t, trainingpredictions = pred) 
	
	return mlmodel

end

#-----------------------------------------------------------------------------------#

function traindynamicsubproblemmodel_mlpass(trainfilelist, traincongestion_flag, trainassignment_flag, setbetas)

	w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, actualobj, oldobj, actualimp, sp_workstations, sp_times, sp_pods, sp_orders = initializedatastructures()

	probleminstances = []
	numsubproblems, weight = Dict(), Dict()

	#Add data for each instance
	for file in trainfilelist #Each file contains one instance and some number of subproblems 

		println(file)

		datafile = string(trainingfolder, file)
		data = load(datafile)

		i = parse(Int64, file[last(findfirst("_run", file))+1 : findfirst(".jld2", file)[1]-1 ])

		push!(probleminstances, i)
		numsubproblems[i] = length(data["subproblemdata"])
		mlpass = parse(Int64, file[last(findfirst("_pass", file))+1 : findfirst("_instance", file)[1]-1])

		#Add instance features
		w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features = addinstance(i, data, featurenames, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums)
		#if mlpass <= 1e-4
		#	push!(probleminstances, i+10000)
		#	w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features = addinstance_nocongestion(i+10000, data, featurenames, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums, featurelookup)
		#	numsubproblems[i+10000] = numsubproblems[i]
		#end

		#Add subproblem structure and objectives
		for k in 1:numsubproblems[i]
			spdata = data["subproblemdata"][k]
			sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp = addsubproblem(i, k, spdata, tstep, targetnumpods, targetnumorders, sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp)
			weight[i,k] = mlpass*passweight + 1
			#if mlpass <= 1e-4
			#	sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp = addsubproblem_nocongestion(i+10000, k, spdata, tstep, targetnumpods, targetnumorders, sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp)
			#	weight[i+10000,k] = mlpass*passweight + 1
			#end
		end
		#if mlpass <= 1e-4
		#	numsubproblems[i] = 2*numsubproblems[i]
		#end

	end

	include("scripts/dynamicsubproblem/mlmodels.jl")
	beta_wt, beta_mp, beta_pw, beta_pwt, pred, shifts = linearregressionmodel_nocongestion(sp_workstations, sp_pods, sp_orders, sp_times, actualobj, probleminstances, numsubproblems, featurelookup, featuresforprediction_wt, featuresforprediction_mp, featuresforprediction_pw, featuresforprediction_pwt, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums, setbetas, traincongestion_flag, trainassignment_flag, weight)
	#printmodelsummary(beta_wt, beta_mp, beta_pw, beta_pwt, pred)

	println("------------------------------------------")
	for f in featuresforprediction_wt
		println(featurenames[f], " = ", value(beta_wt[f]))
	end
	println("------------------------------------------")
	for f in featuresforprediction_mp
		println(featurenames[f], " = ", value(beta_mp[f]))
	end
	println("------------------------------------------")
	for f in featuresforprediction_pw
		println(featurenames[f], " = ", value(beta_pw[f]))
	end
	println("------------------------------------------")
	for f in featuresforprediction_pwt
		println(featurenames[f], " = ", value(beta_pwt[f]))
	end
	println("------------------------------------------")
	println("beta_m = ", shifts[1])
	println("beta_p = ", shifts[2])
	println("beta_w = ", shifts[3])
	println("beta_t = ", shifts[4])
	println("------------------------------------------")

	totalAE, totalSE, tss, validsps = 0, 0, 0, 0
	meanobj = mean(values(actualobj))
	for i in probleminstances, k in 1:numsubproblems[i]
		if actualobj[i,k] > 1e-4
			totalAE += abs(pred[i,k] - actualobj[i,k])
			totalSE += (pred[i,k] - actualobj[i,k])^2
			tss += (meanobj - actualobj[i,k])^2
			validsps += 1
		end
	end

	println("Training mean absolute error = ", totalAE / validsps)
	println("Training mean squared error = ", totalSE / validsps)
	println("Training R_squared = ", 1 - totalSE/tss)

	mlmodel = (features = allfeatures, featurenames = featurenames, beta_wt = beta_wt, beta_mp = beta_mp, beta_pw = beta_pw, beta_pwt = beta_pwt, shifts = shifts, trainingpredictions = pred) 

	transformations = 1:3
	beta_wt2, beta_mp2, beta_pw2, beta_pwt2, pred2, shifts2 = linearregressionmodel_shift_transform(sp_workstations, sp_pods, sp_orders, sp_times, actualobj, probleminstances, numsubproblems, featurelookup, featuresforprediction_wt, featuresforprediction_mp, featuresforprediction_pw, featuresforprediction_pwt, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums, setbetas, traincongestion_flag, trainassignment_flag, weight)
	#printmodelsummary(beta_wt, beta_mp, beta_pw, beta_pwt, pred)

	function printfeature(beta, f)
		if sum(abs(beta[f,tr]) for tr in 2:length(transformations)) > 1e-4
			if abs(beta[f,2]) > 1e-4
				println("log( ", featurenames[f], " )", " = ", value(beta[f,2]))
			elseif abs(beta[f,3]) > 1e-4
				println("( ", featurenames[f]," ) ** 2", " = ", value(beta[f,3]))
			end
		else 
			println(featurenames[f], " = ", value(beta[f,1]))
		end
	end

	println("------------------------------------------")
	for f in featuresforprediction_wt
		printfeature(beta_wt2, f)
	end
	println("------------------------------------------")
	for f in featuresforprediction_mp
		printfeature(beta_mp2, f)
	end
	println("------------------------------------------")
	for f in featuresforprediction_pw
		printfeature(beta_pw2, f)
	end
	println("------------------------------------------")
	for f in featuresforprediction_pwt
		printfeature(beta_pwt2, f)
	end
	println("------------------------------------------")
	println("beta_m = ", shifts2[1])
	println("beta_p = ", shifts2[2])
	println("beta_w = ", shifts2[3])
	println("beta_t = ", shifts2[4])
	println("------------------------------------------")

	totalAE, totalSE, tss, validsps = 0, 0, 0, 0
	meanobj = mean(values(actualobj))
	for i in probleminstances, k in 1:numsubproblems[i]
		if actualobj[i,k] > 1e-4
			totalAE += abs(pred[i,k] - actualobj[i,k])
			totalSE += (pred[i,k] - actualobj[i,k])^2
			tss += (meanobj - actualobj[i,k])^2
			validsps += 1
		end
	end

	println("Training mean absolute error = ", totalAE / validsps)
	println("Training mean squared error = ", totalSE / validsps)
	println("Training R_squared = ", 1 - totalSE/tss)

	mlmodel2 = (features = allfeatures, featurenames = featurenames, beta_wt = beta_wt2, beta_mp = beta_mp2, beta_pw = beta_pw2, beta_pwt = beta_pwt2, shifts = shifts2, trainingpredictions = pred2) 

	return mlmodel, mlmodel2

end

#-----------------------------------------------------------------------------------#

function printmodelsummary(beta_wt, beta_mp, beta_pw, beta_pwt, pred, featurenames, featuresforprediction_wt, featuresforprediction_mp, featuresforprediction_pw, featuresforprediction_pwt)
	
	println("------------------------------------------")
	for f in featuresforprediction_wt
		println(featurenames[f], " = ", value(beta_wt[f]))
	end
	println("------------------------------------------")
	for f in featuresforprediction_mp
		println(featurenames[f], " = ", value(beta_mp[f]))
	end
	println("------------------------------------------")
	for f in featuresforprediction_pw
		println(featurenames[f], " = ", value(beta_pw[f]))
	end
	println("------------------------------------------")
	for f in featuresforprediction_pwt
		println(featurenames[f], " = ", value(beta_pwt[f]))
	end
	println("------------------------------------------")

	totalAE, totalSE, tss, validsps = 0, 0, 0, 0
	meanobj = mean(values(actualobj))
	for i in probleminstances, k in 1:numsubproblems[i]
		if actualobj[i,k] > 1e-4
			totalAE += abs(pred[i,k] - actualobj[i,k])
			totalSE += (pred[i,k] - actualobj[i,k])^2
			tss += (meanobj - actualobj[i,k])^2
			validsps += 1
		end
	end

	println("Training mean absolute error = ", totalAE / validsps)
	println("Training mean squared error = ", totalSE / validsps)
	println("Training R_squared = ", 1 - totalSE/tss)

end

#-----------------------------------------------------------------------------------#

function predictdynamicsubproblemmodel(testfilelist, mlmodel)

	w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, actualobj, oldobj, actualimp, sp_workstations, sp_times, sp_pods, sp_orders = initializedatastructures()
	testinstances = []
	numsubproblems = Dict()
	sp_start, sp_end, sp_ws, sp_podcnt, sp_ordercnt = Dict(), Dict(), Dict(), Dict(), Dict() 
	
	for file in testfilelist

		datafile = string(trainingfolder, file)
		data = load(datafile)

		i = parse(Int64, file[last(findfirst("_run", file))+1 : findfirst(".jld2", file)[1]-1 ])
		push!(testinstances, i)
		numsubproblems[i] = length(data["subproblemdata"])

		orders, workstations, pods, itemson, podswith, horizon, tstep = data["instance_features"]["orders"], data["instance_features"]["workstations"], data["instance_features"]["pods"], data["instance_features"]["itemson"], data["instance_features"]["podswith"], data["instance_features"]["horizon"], data["instance_features"]["tstep"]

		#Add instance features
		w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features = addinstance(i, data, featurenames, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums)

		#Add subproblem structure and objectives
		for k in 1:numsubproblems[i]
			spdata = data["subproblemdata"][k]
			sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp = addsubproblem(i, k, spdata, tstep, targetnumpods, targetnumorders, sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp)
		end

		for k in 1:numsubproblems[i]
			spdata = data["subproblemdata"][k]
			sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp = addsubproblem(i, k, spdata, tstep, targetnumpods, targetnumorders, sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp)

			#Reporting
			#sp_start[i,k] = spdata["sp_tstart"]
			#sp_end[i,k] = spdata["sp_tend"]
			#if length(spdata["sp_workstations"]) == 1
			#	sp_ws[i,k] = spdata["sp_workstations"][1]
			#else
			#	sp_ws[i,k] = "15, 16"
			#end
			#sp_podcnt[i,k] = length(spdata["impact_pods"])
			#sp_ordercnt[i,k] = length(spdata["impact_orders"])
		end

	end

	testpredictions = calcpredictionfortestset_shift(testinstances, mlmodel.beta_wt, mlmodel.beta_mp, mlmodel.beta_pw, mlmodel.beta_pwt, mlmodel.shifts, numsubproblems, sp_workstations, sp_times, sp_pods, sp_orders, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums)

	return testpredictions, actualobj, oldobj, actualimp, testinstances, numsubproblems

end

#-----------------------------------------------------------------------------------#

function predictdynamicsubproblemmodel_compare(testfilelist, mlmodel)

	w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, actualobj, oldobj, actualimp, sp_workstations, sp_times, sp_pods, sp_orders = initializedatastructures()
	testinstances = []
	numsubproblems = Dict()
	sp_start, sp_end, sp_ws, sp_podcnt, sp_ordercnt = Dict(), Dict(), Dict(), Dict(), Dict() 
	
	for file in testfilelist

		datafile = string(trainingfolder, file)
		data = load(datafile)

		i = parse(Int64, file[last(findfirst("_run", file))+1 : findfirst(".jld2", file)[1]-1 ])
		push!(testinstances, i)
		numsubproblems[i] = length(data["subproblemdata"])

		orders, workstations, pods, itemson, podswith, horizon, tstep = data["instance_features"]["orders"], data["instance_features"]["workstations"], data["instance_features"]["pods"], data["instance_features"]["itemson"], data["instance_features"]["podswith"], data["instance_features"]["horizon"], data["instance_features"]["tstep"]

		#Add instance features
		w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features = addinstance(i, data, featurenames, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums)

		#Add subproblem structure and objectives
		for k in 1:numsubproblems[i]
			spdata = data["subproblemdata"][k]
			sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp = addsubproblem(i, k, spdata, tstep, targetnumpods, targetnumorders, sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp)
		end

		for k in 1:numsubproblems[i]
			spdata = data["subproblemdata"][k]
			sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp = addsubproblem(i, k, spdata, tstep, targetnumpods, targetnumorders, sp_workstations, sp_times, sp_pods, sp_orders, actualobj, oldobj, actualimp)

			#Reporting
			#sp_start[i,k] = spdata["sp_tstart"]
			#sp_end[i,k] = spdata["sp_tend"]
			#if length(spdata["sp_workstations"]) == 1
			#	sp_ws[i,k] = spdata["sp_workstations"][1]
			#else
			#	sp_ws[i,k] = "15, 16"
			#end
			#sp_podcnt[i,k] = length(spdata["impact_pods"])
			#sp_ordercnt[i,k] = length(spdata["impact_orders"])
		end

	end

	testpredictions = calcpredictionfortestset_shift(testinstances, mlmodel.beta_wt, mlmodel.beta_mp, mlmodel.beta_pw, mlmodel.beta_pwt, mlmodel.shifts, numsubproblems, sp_workstations, sp_times, sp_pods, sp_orders, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums)

	return testpredictions, actualobj, oldobj, actualimp, testinstances, numsubproblems

end


#-----------------------------------------------------------------------------------#

function printmodelmetrics(pred, actual, probleminstances, actualimp, oldobj, numsubproblems)
	
	totalAE, totalSE, tss, validsps = 0, 0, 0, 0
	meanobj = mean(values(actual))
	splist = []
	for i in probleminstances, k in 1:numsubproblems[i]
		push!(splist, (i,k))
		if actual[i,k] > 1e-4
			totalAE += abs(pred[i,k] - actual[i,k])
			totalSE += (pred[i,k] - actual[i,k])^2
			tss += (meanobj - actual[i,k])^2
			validsps += 1
		end
	end 

	println("Testing mean absolute error = ", totalAE / validsps)
	println("Testing mean squared error = ", totalSE / validsps)
	println("Testing R_squared = ", 1 - totalSE/tss)

	numsubproblems = length(splist)

	sortedimp = sort(collect(actualimp), by = x->x[2])
	top10_act = [(i,k) for (i,k) in splist if actualimp[i,k] >= sortedimp[numsubproblems - convert(Int, ceil(numsubproblems*0.1))][2]]
	top5_act = [(i,k) for (i,k) in splist if actualimp[i,k] >= sortedimp[numsubproblems - convert(Int, ceil(numsubproblems*0.05))][2]]
	top1_act = [(i,k) for (i,k) in splist if actualimp[i,k] >= sortedimp[numsubproblems - convert(Int, ceil(numsubproblems*0.01))][2]]
	top01_act = [(i,k) for (i,k) in splist if actualimp[i,k] >= sortedimp[numsubproblems - convert(Int, ceil(numsubproblems*0.001))][2]]

	predimp = [pred[sp] - oldobj[sp] for sp in splist]
	top10_pred = [(i,k) for (i,k) in splist if pred[i,k] - oldobj[i,k] >= sort(predimp)[numsubproblems - convert(Int, ceil(numsubproblems*0.1))]]
	top5_pred = [(i,k) for (i,k) in splist if pred[i,k] - oldobj[i,k] >= sort(predimp)[numsubproblems - convert(Int, ceil(numsubproblems*0.05))]]
	top1_pred = [(i,k) for (i,k) in splist if pred[i,k] - oldobj[i,k] >= sort(predimp)[numsubproblems - convert(Int, ceil(numsubproblems*0.01))]]
	top01_pred = [(i,k) for (i,k) in splist if pred[i,k] - oldobj[i,k] >= sort(predimp)[numsubproblems - convert(Int, ceil(numsubproblems*0.001))]]

	println("Top 10 pct identification accuracy = ", length(intersect(top10_act, top10_pred)) / length(top10_pred) )
	println("Top 5 pct identification accuracy = ", length(intersect(top5_act, top5_pred)) / length(top5_pred) )
	println("Top 1 pct identification accuracy = ", length(intersect(top1_act, top1_pred)) / length(top1_pred) )
	println("Top 0.1 pct identification accuracy = ", length(intersect(top01_act, top01_pred)) / length(top01_pred) )

end

