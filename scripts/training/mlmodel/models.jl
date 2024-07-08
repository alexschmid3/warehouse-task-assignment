

function linearregressionmodel(sp_workstations, sp_pods, sp_orders, sp_times, actualobj, probleminstances, numsubproblems, featurelookup, featuresforprediction_wt, featuresforprediction_mp, featuresforprediction_pw, featuresforprediction_pwt, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums, setbetas, traincongestion_flag, trainassignment_flag, weight)

	#consistency_flag, x, y, z, v, q, actualobj = 1, length(workstationgrouplist), x_k, y_k, z_k, v_k, q_k, obj_k

	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", 60*60) # iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", 0)
	set_optimizer_attribute(model, "IntegralityFocus", 1)

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
	@objective(model, Min, sum(sum(weight[i,k] * err[i,k] for k in 1:numsubproblems[i]) for i in probleminstances))

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

	#Prediction error
	#@constraint(model, prederror1[i = probleminstances, k = 1:numsubproblems[i]], err[i,k] >= pred[i,k] - actualobj[i,k])
	#@constraint(model, prederror2[i = probleminstances, k = 1:numsubproblems[i]], err[i,k] >= actualobj[i,k] - pred[i,k])
	for i in probleminstances, k in 1:numsubproblems[i]
		if actualobj[i,k] > 1e-4
			@constraint(model, err[i,k] >= pred[i,k] - actualobj[i,k])
			@constraint(model, err[i,k] >= actualobj[i,k] - pred[i,k])
		else
			@constraint(model, err[i,k] == 0)
		end
	end

	#Consistency / interpretability of predictions
	#@constraint(model, otherinstancesareworse[i = probleminstances, k = 1:numsubproblems[i]], pred[i, k] <= pred[i, k_star])
	#@constraint(model, nonzerofeatures_wt_1[f in featuresforprediction_wt], beta_wt[f] <= 10000 * nonzero_wt[f])
	#@constraint(model, nonzerofeatures_wt_2[f in featuresforprediction_wt], beta_wt[f] >= -10000 * nonzero_wt[f])
	@constraint(model, nonzerofeatures_mp_1[f in featuresforprediction_mp], beta_mp[f] <= 10000 * nonzero_mp[f])
	@constraint(model, nonzerofeatures_mp_2[f in featuresforprediction_mp], beta_mp[f] >= -10000 * nonzero_mp[f])
	@constraint(model, nonzerofeatures_pw_1[f in featuresforprediction_pw], beta_pw[f] <= 10000 * nonzero_pw[f])
	@constraint(model, nonzerofeatures_pw_2[f in featuresforprediction_pw], beta_pw[f] >= -10000 * nonzero_pw[f])
	#@constraint(model, nonzerofeatures_pwt_1[f in featuresforprediction_pwt], beta_pwt[f] <= 10000 * nonzero_pwt[f])
	#@constraint(model, nonzerofeatures_pwt_2[f in featuresforprediction_pwt], beta_pwt[f] >= -10000 * nonzero_pwt[f])
	@constraint(model, nonzero_mp[featurelookup["mp_itemoverlap"]] + nonzero_mp[featurelookup["mp_itemoverlappct"]] <= 1)
	@constraint(model, nonzero_mp[featurelookup["m_ordersize"]] + nonzero_mp[featurelookup["m_oneitemflag"]] <= 1)
	@constraint(model, nonzero_mp[featurelookup["mp_overlapitemavginv"]] + nonzero_mp[featurelookup["mp_overlapitemavgaltpods"]] <= 1)
	@constraint(model, nonzero_pw[featurelookup["pw_distance"]] + nonzero_pw[featurelookup["pw_xsteps"]] <= 1)
	@constraint(model, nonzero_pw[featurelookup["pw_distance"]] + nonzero_pw[featurelookup["pw_ysteps"]] <= 1)

	@constraint(model, beta_mp[featurelookup["p_pctusefulitems"]] == 0)
	@constraint(model, beta_pw[featurelookup["p_pctusefulitems"]] == 0)
	@constraint(model, beta_pwt[featurelookup["p_avgdisttostations"]] == 0)
	@constraint(model, beta_pw[featurelookup["p_avgdisttostations"]] == 0)
	@constraint(model, beta_pw[featurelookup["w_centrality"]] == 0)

	@constraint(model, beta_pwt[featurelookup["t_horizonrelativetime"]] == 0)
	@constraint(model, beta_wt[featurelookup["t_horizonrelativetime"]] == 0)

	#Isolate training
	setbeta_wt, setbeta_mp, setbeta_pw, setbeta_pwt = setbetas
	if trainassignment_flag == 1 
		@constraint(model, force_wt[f = congestionfeatures_wt], beta_wt[f] == 0)
		@constraint(model, force_pwt[f = congestionfeatures_pwt], beta_pwt[f] == 0)
	elseif traincongestion_flag == 1
		@constraint(model, force_wt[f = setdiff(featuresforprediction_wt, congestionfeatures_wt)], beta_wt[f] == setbeta_wt[f])
		@constraint(model, force_mp[f = featuresforprediction_mp], beta_mp[f] == setbeta_mp[f])
		@constraint(model, force_pw[f = featuresforprediction_pw], beta_pw[f] == setbeta_pw[f])
		@constraint(model, force_pwt[f = setdiff(featuresforprediction_pwt, congestionfeatures_pwt)], beta_pwt[f] == setbeta_pwt[f])
	end

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

function linearregressionmodel_shift(sp_workstations, sp_pods, sp_orders, sp_times, actualobj, probleminstances, numsubproblems, featurelookup, featuresforprediction_wt, featuresforprediction_mp, featuresforprediction_pw, featuresforprediction_pwt, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums, setbetas, traincongestion_flag, trainassignment_flag, weight)

	#consistency_flag, x, y, z, v, q, actualobj = 1, length(workstationgrouplist), x_k, y_k, z_k, v_k, q_k, obj_k

	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", 60*60) # iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", 1)
	set_optimizer_attribute(model, "IntegralityFocus", 1)

	#Variables
	@variable(model, pred[i = probleminstances, k = 1:numsubproblems[i]])
	@variable(model, err[i = probleminstances, k = 1:numsubproblems[i]])
	@variable(model, beta_wt[f = featuresforprediction_wt])
	@variable(model, beta_mp[f = featuresforprediction_mp])
	@variable(model, beta_pw[f = featuresforprediction_pw])
	@variable(model, beta_pwt[f = featuresforprediction_pwt])

	@variable(model, beta_t)
	@variable(model, beta_p)
	@variable(model, beta_m)
	@variable(model, beta_w)

	@variable(model, nonzero_wt[f = featuresforprediction_wt], Bin)
	@variable(model, nonzero_mp[f = featuresforprediction_mp], Bin)
	@variable(model, nonzero_pw[f = featuresforprediction_pw], Bin)
	@variable(model, nonzero_pwt[f = featuresforprediction_pwt], Bin)
	@variable(model, chooseone, Bin)

	#Objective
	@objective(model, Min, sum(sum(weight[i,k] * err[i,k] for k in 1:numsubproblems[i]) for i in probleminstances))

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
		#Bases
		+ sum(beta_t for t in sp_times[i,k]) 
		+ sum(beta_p for t in sp_pods[i,k])
		+ sum(beta_m for t in sp_orders[i,k])
		+ sum(beta_w for t in sp_workstations[i,k])
		)

	#Prediction error
	for i in probleminstances, k in 1:numsubproblems[i]
		if actualobj[i,k] > 1e-4
			@constraint(model, err[i,k] >= pred[i,k] - actualobj[i,k])
			@constraint(model, err[i,k] >= actualobj[i,k] - pred[i,k])
		else
			@constraint(model, err[i,k] == 0)
		end
	end

	#Consistency / interpretability of predictions                                                                    
	@constraint(model, nonzerofeatures_wt_1[f in featuresforprediction_wt], beta_wt[f] <= 10000 * nonzero_wt[f])
	@constraint(model, nonzerofeatures_wt_2[f in featuresforprediction_wt], beta_wt[f] >= -10000 * nonzero_wt[f])
	@constraint(model, nonzerofeatures_mp_1[f in featuresforprediction_mp], beta_mp[f] <= 10000 * nonzero_mp[f])
	@constraint(model, nonzerofeatures_mp_2[f in featuresforprediction_mp], beta_mp[f] >= -10000 * nonzero_mp[f])
	@constraint(model, nonzerofeatures_pw_1[f in featuresforprediction_pw], beta_pw[f] <= 10000 * nonzero_pw[f])
	@constraint(model, nonzerofeatures_pw_2[f in featuresforprediction_pw], beta_pw[f] >= -10000 * nonzero_pw[f])
	@constraint(model, nonzerofeatures_pwt_1[f in featuresforprediction_pwt], beta_pwt[f] <= 10000 * nonzero_pwt[f])
	@constraint(model, nonzerofeatures_pwt_2[f in featuresforprediction_pwt], beta_pwt[f] >= -10000 * nonzero_pwt[f])
	@constraint(model, nonzero_mp[featurelookup["mp_itemoverlap"]] + nonzero_mp[featurelookup["mp_itemoverlappct"]] <= 1)
	@constraint(model, nonzero_mp[featurelookup["m_ordersize"]] + nonzero_mp[featurelookup["m_oneitemflag"]] <= 1)
	@constraint(model, nonzero_mp[featurelookup["mp_overlapitemavginv"]] + nonzero_mp[featurelookup["mp_overlapitemavgaltpods"]] <= 1)
	@constraint(model, nonzero_pw[featurelookup["pw_distance"]] + nonzero_pw[featurelookup["pw_xsteps"]] <= 1)
	@constraint(model, nonzero_pw[featurelookup["pw_distance"]] + nonzero_pw[featurelookup["pw_ysteps"]] <= 1)
	@constraint(model, nonzero_wt[featurelookup["wt_avglocalcong"]] + nonzero_wt[featurelookup["wt_maxlocalcong"]] <= 1)
	@constraint(model, nonzero_wt[featurelookup["wt_avghyperlocalcong"]] + nonzero_wt[featurelookup["wt_maxhyperlocalcong"]] <= 1)
	@constraint(model, nonzero_wt[featurelookup["w_centrality"]] + nonzero_pwt[featurelookup["w_centrality"]] <= 1)

	@constraint(model, beta_mp[featurelookup["p_pctusefulitems"]] == 0)
	@constraint(model, beta_pw[featurelookup["p_pctusefulitems"]] == 0)
	@constraint(model, beta_pwt[featurelookup["p_avgdisttostations"]] == 0)
	@constraint(model, beta_pw[featurelookup["p_avgdisttostations"]] == 0)
	@constraint(model, beta_pw[featurelookup["w_centrality"]] == 0)

	@constraint(model, beta_pwt[featurelookup["t_horizonrelativetime"]] == 0)
	@constraint(model, beta_wt[featurelookup["t_horizonrelativetime"]] == 0)
	@constraint(model, beta_wt[featurelookup["t_timeofday"]] == 0)
	@constraint(model, beta_pwt[featurelookup["t_timeofday"]] == 0)
	@constraint(model, beta_pwt[featurelookup["pwt_existingcong"]] == 0.0)
	
	@constraint(model, beta_m == 0)
	@constraint(model, beta_p == 0)

	@constraint(model, beta_wt[featurelookup["wt_queuecong"]] == 0.0)
	@constraint(model, beta_wt[featurelookup["wt_avglocalcong"]] ==  0.0)
	@constraint(model, beta_wt[featurelookup["wt_maxlocalcong"]] == 0.0)
	@constraint(model, beta_wt[featurelookup["wt_avghyperlocalcong"]] == 0.0)
	@constraint(model, beta_wt[featurelookup["wt_maxhyperlocalcong"]] == 0.0)
	
	#@constraint(model, beta_wt[featurelookup["wt_maxlocalcong"]] == 0)
	#@constraint(model, beta_wt[featurelookup["wt_queuecong"]] == 0)
	#@constraint(model, beta_wt[featurelookup["wt_avglocalcong"]] == 0)
	#@constraint(model, beta_wt[featurelookup["wt_avghyperlocalcong"]] == 0)
	#@constraint(model, beta_wt[featurelookup["wt_maxhyperlocalcong"]] == 0)
	#@constraint(model, beta_pwt[featurelookup["pwt_existingcong"]] == 0)

	#Isolate training
	setbeta_wt, setbeta_mp, setbeta_pw, setbeta_pwt = setbetas
	if trainassignment_flag == 1 
		@constraint(model, force_wt[f = congestionfeatures_wt], beta_wt[f] == 0)
		@constraint(model, force_pwt[f = congestionfeatures_pwt], beta_pwt[f] == 0)
	elseif traincongestion_flag == 1
		@constraint(model, force_wt[f = setdiff(featuresforprediction_wt, congestionfeatures_wt)], beta_wt[f] == setbeta_wt[f])
		@constraint(model, force_mp[f = featuresforprediction_mp], beta_mp[f] == setbeta_mp[f])
		@constraint(model, force_pw[f = featuresforprediction_pw], beta_pw[f] == setbeta_pw[f])
		@constraint(model, force_pwt[f = setdiff(featuresforprediction_pwt, congestionfeatures_pwt)], beta_pwt[f] == setbeta_pwt[f])
	end

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

	return value.(beta_wt), value.(beta_mp), value.(beta_pw), value.(beta_pwt), value.(pred), [value(beta_m), value(beta_p), value(beta_w), value(beta_t)]

end

#-----------------------------------------------------------------------------------------------#

function transformfeat(tr, num)
	if tr == 1
		return num
	elseif (tr == 2) * (num > 0)
		return log(num)
	elseif  (tr == 2) 
		return num
	elseif tr == 3
		return num^2
	end
end

#-----------------------------------------------------------------------------------------------#

function linearregressionmodel_shift_transform(sp_workstations, sp_pods, sp_orders, sp_times, actualobj, probleminstances, numsubproblems, featurelookup, featuresforprediction_wt, featuresforprediction_mp, featuresforprediction_pw, featuresforprediction_pwt, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums, setbetas, traincongestion_flag, trainassignment_flag, weight)

	#consistency_flag, x, y, z, v, q, actualobj = 1, length(workstationgrouplist), x_k, y_k, z_k, v_k, q_k, obj_k

	transformations = 1:3 #1-none, 2-log, 3-squared
		
	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", 60*60) # iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", 1)
	set_optimizer_attribute(model, "IntegralityFocus", 1)

	#Variables
	@variable(model, pred[i = probleminstances, k = 1:numsubproblems[i]])
	@variable(model, err[i = probleminstances, k = 1:numsubproblems[i]])

	@variable(model, beta_wt[f = featuresforprediction_wt, t = transformations])
	@variable(model, beta_mp[f = featuresforprediction_mp, t = transformations])
	@variable(model, beta_pw[f = featuresforprediction_pw, t = transformations])
	@variable(model, beta_pwt[f = featuresforprediction_pwt, t = transformations])

	@variable(model, beta_t)
	@variable(model, beta_p)
	@variable(model, beta_m)
	@variable(model, beta_w)

	@variable(model, nonzero_wt[f = featuresforprediction_wt, tr = transformations], Bin)
	@variable(model, nonzero_mp[f = featuresforprediction_mp, tr = transformations], Bin)
	@variable(model, nonzero_pw[f = featuresforprediction_pw, tr = transformations], Bin)
	@variable(model, nonzero_pwt[f = featuresforprediction_pwt, tr = transformations], Bin)
	@variable(model, chooseone, Bin)

	#Regularization
	@variable(model, abs_beta_wt[f = featuresforprediction_wt, tr = transformations])
	@variable(model, abs_beta_mp[f = featuresforprediction_mp, tr = transformations])
	@variable(model, abs_beta_pw[f = featuresforprediction_pw, tr = transformations])
	@variable(model, abs_beta_pwt[f = featuresforprediction_pwt, tr = transformations])

	#Objective
	@objective(model, Min, sum(sum(weight[i,k] * err[i,k] for k in 1:numsubproblems[i]) for i in probleminstances) 
		+ robustnessparam * (
			sum(abs_beta_wt[f,tr] for f in featuresforprediction_wt, tr in transformations) 
			+ sum(abs_beta_mp[f,tr] for f in featuresforprediction_mp, tr in transformations)  
			+ sum(abs_beta_pw[f,tr] for f in featuresforprediction_pw, tr in transformations) 
			+ sum(abs_beta_pwt[f,tr] for f in featuresforprediction_pwt, tr in transformations) 
		))

	#Prediction calculation
	@constraint(model, predictioncalc[i = probleminstances, k = 1:numsubproblems[i]], pred[i,k] == 
		sum(#Workstation-time
			sum(sum(sum(beta_wt[f,tr] * transformfeat(tr, wt_features[i,f,w,t]) for f in wt_featnums) for w in sp_workstations[i,k]) for t in sp_times[i,k])
			+ sum(sum(beta_wt[f,tr] * transformfeat(tr, w_features[i,f,w]) for f in w_featnums) for w in sp_workstations[i,k]) 
			+ sum(sum(beta_wt[f,tr] * transformfeat(tr, t_features[i,f,t]) for f in t_featnums) for t in sp_times[i,k]) 
			#Order-pod
			+ sum(sum(sum(beta_mp[f,tr] * transformfeat(tr, mp_features[i,f,m,p]) for f in mp_featnums) for m in sp_orders[i,k]) for p in sp_pods[i,k]) 
			+ sum(sum(beta_mp[f,tr] * transformfeat(tr, m_features[i,f,m]) for f in m_featnums) for m in sp_orders[i,k]) 
			+ sum(sum(beta_mp[f,tr] * transformfeat(tr, p_features[i,f,p]) for f in p_featnums) for p in sp_pods[i,k]) 
			#Pod-workstation
			+ sum(sum(sum(beta_pw[f,tr] * transformfeat(tr, pw_features[i,f,p,w]) for f in pw_featnums) for w in sp_workstations[i,k]) for p in sp_pods[i,k]) 
			+ sum(sum(beta_pw[f,tr] * transformfeat(tr, p_features[i,f,p]) for f in p_featnums) for p in sp_pods[i,k]) 
			+ sum(sum(beta_pw[f,tr] * transformfeat(tr, w_features[i,f,w]) for f in w_featnums) for w in sp_workstations[i,k]) 
			#Pod-workstation-time
			+ sum(sum(sum(sum(beta_pwt[f,tr] * transformfeat(tr, pwt_features[i,f,p,w,t]) for f in pwt_featnums) for w in sp_workstations[i,k]) for p in sp_pods[i,k]) for t in sp_times[i,k])  
			+ sum(sum(beta_pwt[f,tr] * transformfeat(tr, p_features[i,f,p]) for f in p_featnums) for p in sp_pods[i,k])
			+ sum(sum(beta_pwt[f,tr] * transformfeat(tr, w_features[i,f,w]) for f in w_featnums) for w in sp_workstations[i,k]) 
			+ sum(sum(beta_pwt[f,tr] * transformfeat(tr, t_features[i,f,t]) for f in t_featnums) for t in sp_times[i,k])
		for tr in transformations)
		#Bases
		+ sum(beta_t for t in sp_times[i,k]) 
		+ sum(beta_p for t in sp_pods[i,k])
		+ sum(beta_m for t in sp_orders[i,k])
		+ sum(beta_w for t in sp_workstations[i,k])
		)

	#Prediction error
	for i in probleminstances, k in 1:numsubproblems[i]
		if actualobj[i,k] > 1e-4
			@constraint(model, err[i,k] >= pred[i,k] - actualobj[i,k])
			@constraint(model, err[i,k] >= actualobj[i,k] - pred[i,k])
		else
			@constraint(model, err[i,k] == 0)
		end
	end

	#Regularization
	@constraint(model, absval_beta_wt1[f = featuresforprediction_wt, tr = transformations], abs_beta_wt[f,tr] >= beta_wt[f,tr])
	@constraint(model, absval_beta_wt2[f = featuresforprediction_wt, tr = transformations], abs_beta_wt[f,tr] >= - beta_wt[f,tr])
	@constraint(model, absval_beta_mp1[f = featuresforprediction_mp, tr = transformations], abs_beta_mp[f,tr] >= beta_mp[f,tr])
	@constraint(model, absval_beta_mp2[f = featuresforprediction_mp, tr = transformations], abs_beta_mp[f,tr] >= - beta_mp[f,tr])
	@constraint(model, absval_beta_pw1[f = featuresforprediction_pw, tr = transformations], abs_beta_pw[f,tr] >= beta_pw[f,tr])
	@constraint(model, absval_beta_pw2[f = featuresforprediction_pw, tr = transformations], abs_beta_pw[f,tr] >= - beta_pw[f,tr])
	@constraint(model, absval_beta_pwt1[f = featuresforprediction_pwt, tr = transformations], abs_beta_pwt[f,tr] >= beta_pwt[f,tr])
	@constraint(model, absval_beta_pwt2[f = featuresforprediction_pwt, tr = transformations], abs_beta_pwt[f,tr] >= - beta_pwt[f,tr])

	#Consistency / interpretability of predictions                                                                    
	@constraint(model, nonzerofeatures_wt_1[f in featuresforprediction_wt, tr in transformations], beta_wt[f, tr] <= 10000 * nonzero_wt[f, tr])
	@constraint(model, nonzerofeatures_wt_2[f in featuresforprediction_wt, tr in transformations], beta_wt[f, tr] >= -10000 * nonzero_wt[f, tr])
	@constraint(model, nonzerofeatures_mp_1[f in featuresforprediction_mp, tr in transformations], beta_mp[f, tr] <= 10000 * nonzero_mp[f, tr])
	@constraint(model, nonzerofeatures_mp_2[f in featuresforprediction_mp, tr in transformations], beta_mp[f, tr] >= -10000 * nonzero_mp[f, tr])
	@constraint(model, nonzerofeatures_pw_1[f in featuresforprediction_pw, tr in transformations], beta_pw[f, tr] <= 10000 * nonzero_pw[f, tr])
	@constraint(model, nonzerofeatures_pw_2[f in featuresforprediction_pw, tr in transformations], beta_pw[f, tr] >= -10000 * nonzero_pw[f, tr])
	@constraint(model, nonzerofeatures_pwt_1[f in featuresforprediction_pwt, tr in transformations], beta_pwt[f, tr] <= 10000 * nonzero_pwt[f, tr])
	@constraint(model, nonzerofeatures_pwt_2[f in featuresforprediction_pwt, tr in transformations], beta_pwt[f, tr] >= -10000 * nonzero_pwt[f, tr])
	
	@constraint(model, sum(nonzero_mp[featurelookup["mp_itemoverlap"],tr] + nonzero_mp[featurelookup["mp_itemoverlappct"],tr] for tr in transformations) <= 1)
	@constraint(model, sum(nonzero_mp[featurelookup["m_ordersize"],tr] + nonzero_mp[featurelookup["m_oneitemflag"],tr] for tr in transformations)  <= 1)
	@constraint(model, sum(nonzero_mp[featurelookup["mp_overlapitemavginv"],tr] + nonzero_mp[featurelookup["mp_overlapitemavgaltpods"],tr] for tr in transformations)  <= 1)
	@constraint(model, sum(nonzero_pw[featurelookup["pw_distance"],tr] + nonzero_pw[featurelookup["pw_xsteps"],tr] for tr in transformations)  <= 1)
	@constraint(model, sum(nonzero_pw[featurelookup["pw_distance"],tr] + nonzero_pw[featurelookup["pw_ysteps"],tr] for tr in transformations)  <= 1)
	@constraint(model, sum(nonzero_wt[featurelookup["wt_avglocalcong"],tr] + nonzero_wt[featurelookup["wt_maxlocalcong"],tr] for tr in transformations)  <= 1)
	@constraint(model, sum(nonzero_wt[featurelookup["wt_avghyperlocalcong"],tr] + nonzero_wt[featurelookup["wt_maxhyperlocalcong"],tr] for tr in transformations)  <= 1)
	@constraint(model, sum(nonzero_wt[featurelookup["w_centrality"],tr] + nonzero_pwt[featurelookup["w_centrality"],tr] for tr in transformations)  <= 1)

	#Un-used features
	@constraint(model, [tr in transformations], beta_mp[featurelookup["p_pctusefulitems"],tr] == 0)
	@constraint(model, [tr in transformations], beta_pw[featurelookup["p_pctusefulitems"],tr] == 0)
	@constraint(model, [tr in transformations], beta_pwt[featurelookup["p_avgdisttostations"],tr] == 0)
	@constraint(model, [tr in transformations], beta_pw[featurelookup["p_avgdisttostations"],tr] == 0)
	@constraint(model, [tr in transformations], beta_pw[featurelookup["w_centrality"],tr] == 0)

	@constraint(model, [tr in transformations], beta_pwt[featurelookup["t_horizonrelativetime"],tr] == 0)
	@constraint(model, [tr in transformations], beta_wt[featurelookup["t_horizonrelativetime"],tr] == 0)
	@constraint(model, [tr in transformations], beta_wt[featurelookup["t_timeofday"],tr] == 0)
	@constraint(model, [tr in transformations], beta_pwt[featurelookup["t_timeofday"],tr] == 0)
	@constraint(model, [tr in transformations], beta_pwt[featurelookup["pwt_existingcong"],tr] == 0.0)
	
	@constraint(model, beta_m == 0)
	@constraint(model, beta_p == 0)

	@constraint(model, [tr in transformations], beta_wt[featurelookup["wt_queuecong"],tr] == 0.0)
	@constraint(model, [tr in transformations], beta_wt[featurelookup["wt_avglocalcong"],tr] ==  0.0)
	@constraint(model, [tr in transformations], beta_wt[featurelookup["wt_maxlocalcong"],tr] == 0.0)
	@constraint(model, [tr in transformations], beta_wt[featurelookup["wt_avghyperlocalcong"],tr] == 0.0)
	@constraint(model, [tr in transformations], beta_wt[featurelookup["wt_maxhyperlocalcong"],tr] == 0.0)
	
	#@constraint(model, [tr in transformations], beta_wt[featurelookup["wt_maxlocalcong"],tr] == 0)
	#@constraint(model, [tr in transformations], beta_wt[featurelookup["wt_queuecong"],tr] == 0)
	#@constraint(model, [tr in transformations], beta_wt[featurelookup["wt_avglocalcong"],tr] == 0)
	#@constraint(model, [tr in transformations], beta_wt[featurelookup["wt_avghyperlocalcong"],tr] == 0)
	#@constraint(model, [tr in transformations], beta_wt[featurelookup["wt_maxhyperlocalcong"],tr] == 0)
	#@constraint(model, [tr in transformations], beta_pwt[featurelookup["pwt_existingcong"], tr] == 0)

	#Non-linear trasnformations
	@constraint(model, nonlinear_wt[f = featuresforprediction_wt], sum(nonzero_wt[f,tr] for tr in transformations) <= 1)
	@constraint(model, nonlinear_mp[f = featuresforprediction_mp], sum(nonzero_mp[f,tr] for tr in transformations) <= 1)
	@constraint(model, nonlinear_pw[f = featuresforprediction_pw], sum(nonzero_pw[f,tr] for tr in transformations) <= 1)
	@constraint(model, nonlinear_pwt[f = featuresforprediction_pwt], sum(nonzero_pwt[f,tr] for tr in transformations) <= 1)

	#Isolate training
	setbeta_wt, setbeta_mp, setbeta_pw, setbeta_pwt = setbetas
	if trainassignment_flag == 1 
		@constraint(model, force_wt[f = congestionfeatures_wt], beta_wt[f] == 0)
		@constraint(model, force_pwt[f = congestionfeatures_pwt], beta_pwt[f] == 0)
	elseif traincongestion_flag == 1
		@constraint(model, force_wt[f = setdiff(featuresforprediction_wt, congestionfeatures_wt)], beta_wt[f] == setbeta_wt[f])
		@constraint(model, force_mp[f = featuresforprediction_mp], beta_mp[f] == setbeta_mp[f])
		@constraint(model, force_pw[f = featuresforprediction_pw], beta_pw[f] == setbeta_pw[f])
		@constraint(model, force_pwt[f = setdiff(featuresforprediction_pwt, congestionfeatures_pwt)], beta_pwt[f] == setbeta_pwt[f])
	end

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

	return value.(beta_wt), value.(beta_mp), value.(beta_pw), value.(beta_pwt), value.(pred), [value(beta_m), value(beta_p), value(beta_w), value(beta_t)]

end

#-----------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------#

function linearregressionmodel_wt(sp_workstations, sp_pods, sp_orders, sp_times, actualobj, probleminstances, numsubproblems, featurelookup, featuresforprediction_wt, featuresforprediction_mp, featuresforprediction_pw, featuresforprediction_pwt, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums, setbetas, traincongestion_flag, trainassignment_flag, weight)

	#consistency_flag, x, y, z, v, q, actualobj = 1, length(workstationgrouplist), x_k, y_k, z_k, v_k, q_k, obj_k

	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", 60*60) # iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", 1)
	set_optimizer_attribute(model, "IntegralityFocus", 1)

	#Variables
	@variable(model, pred[i = probleminstances, k = 1:numsubproblems[i]])
	@variable(model, err[i = probleminstances, k = 1:numsubproblems[i]])
	@variable(model, beta_wt[f = featuresforprediction_wt])
	@variable(model, beta_mp[f = featuresforprediction_mp])
	@variable(model, beta_pw[f = featuresforprediction_pw])
	@variable(model, beta_pwt[f = featuresforprediction_pwt])

	@variable(model, beta_t)
	@variable(model, beta_p)
	@variable(model, beta_m)
	@variable(model, beta_w)

	@variable(model, nonzero_wt[f = featuresforprediction_wt], Bin)
	@variable(model, nonzero_mp[f = featuresforprediction_mp], Bin)
	@variable(model, nonzero_pw[f = featuresforprediction_pw], Bin)
	@variable(model, nonzero_pwt[f = featuresforprediction_pwt], Bin)
	@variable(model, chooseone, Bin)

	#Objective
	@objective(model, Min, sum(sum(weight[i,k] * err[i,k] for k in 1:numsubproblems[i]) for i in probleminstances))

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
		#Bases
		+ sum(beta_t for t in sp_times[i,k]) 
		+ sum(beta_p for t in sp_pods[i,k])
		+ sum(beta_m for t in sp_orders[i,k])
		+ sum(beta_w for t in sp_workstations[i,k])
		)

	#Prediction error
	#@constraint(model, prederror1[i = probleminstances, k = 1:numsubproblems[i]], err[i,k] >= pred[i,k] - actualobj[i,k])
	#@constraint(model, prederror2[i = probleminstances, k = 1:numsubproblems[i]], err[i,k] >= actualobj[i,k] - pred[i,k])
	for i in probleminstances, k in 1:numsubproblems[i]
		if actualobj[i,k] > 1e-4
			@constraint(model, err[i,k] >= pred[i,k] - actualobj[i,k])
			@constraint(model, err[i,k] >= actualobj[i,k] - pred[i,k])
		else
			@constraint(model, err[i,k] == 0)
		end
	end

	#Consistency / interpretability of predictions
	#@constraint(model, otherinstancesareworse[i = probleminstances, k = 1:numsubproblems[i]], pred[i, k] <= pred[i, k_star])
	@constraint(model, nonzerofeatures_wt_1[f in featuresforprediction_wt], beta_wt[f] <= 10000 * nonzero_wt[f])
	@constraint(model, nonzerofeatures_wt_2[f in featuresforprediction_wt], beta_wt[f] >= -10000 * nonzero_wt[f])
	@constraint(model, nonzerofeatures_mp_1[f in featuresforprediction_mp], beta_mp[f] <= 10000 * nonzero_mp[f])
	@constraint(model, nonzerofeatures_mp_2[f in featuresforprediction_mp], beta_mp[f] >= -10000 * nonzero_mp[f])
	@constraint(model, nonzerofeatures_pw_1[f in featuresforprediction_pw], beta_pw[f] <= 10000 * nonzero_pw[f])
	@constraint(model, nonzerofeatures_pw_2[f in featuresforprediction_pw], beta_pw[f] >= -10000 * nonzero_pw[f])
	@constraint(model, nonzerofeatures_pwt_1[f in featuresforprediction_pwt], beta_pwt[f] <= 10000 * nonzero_pwt[f])
	@constraint(model, nonzerofeatures_pwt_2[f in featuresforprediction_pwt], beta_pwt[f] >= -10000 * nonzero_pwt[f])
	@constraint(model, nonzero_mp[featurelookup["mp_itemoverlap"]] + nonzero_mp[featurelookup["mp_itemoverlappct"]] <= 1)
	@constraint(model, nonzero_mp[featurelookup["m_ordersize"]] + nonzero_mp[featurelookup["m_oneitemflag"]] <= 1)
	@constraint(model, nonzero_mp[featurelookup["mp_overlapitemavginv"]] + nonzero_mp[featurelookup["mp_overlapitemavgaltpods"]] <= 1)
	@constraint(model, nonzero_pw[featurelookup["pw_distance"]] + nonzero_pw[featurelookup["pw_xsteps"]] <= 1)
	@constraint(model, nonzero_pw[featurelookup["pw_distance"]] + nonzero_pw[featurelookup["pw_ysteps"]] <= 1)
	@constraint(model, nonzero_wt[featurelookup["wt_avglocalcong"]] + nonzero_wt[featurelookup["wt_maxlocalcong"]] <= 1)
	@constraint(model, nonzero_wt[featurelookup["wt_avghyperlocalcong"]] + nonzero_wt[featurelookup["wt_maxhyperlocalcong"]] <= 1)
	@constraint(model, nonzero_wt[featurelookup["w_centrality"]] + nonzero_pwt[featurelookup["w_centrality"]] <= 1)

	@constraint(model, beta_mp[featurelookup["p_pctusefulitems"]] == 0)
	@constraint(model, beta_pw[featurelookup["p_pctusefulitems"]] == 0)
	@constraint(model, beta_pwt[featurelookup["p_avgdisttostations"]] == 0)
	@constraint(model, beta_pw[featurelookup["p_avgdisttostations"]] == 0)
	@constraint(model, beta_pw[featurelookup["w_centrality"]] == 0)

	@constraint(model, beta_pwt[featurelookup["t_horizonrelativetime"]] == 0)
	@constraint(model, beta_wt[featurelookup["t_horizonrelativetime"]] == 0)
	@constraint(model, beta_wt[featurelookup["t_timeofday"]] == 0)
	@constraint(model, beta_pwt[featurelookup["t_timeofday"]] == 0)
	@constraint(model, beta_pwt[featurelookup["pwt_existingcong"]] == 0.0)
	
	@constraint(model, beta_m == 0)
	@constraint(model, beta_p == 0)
	
	#@constraint(model, beta_wt[featurelookup["wt_maxlocalcong"]] == 0)
	#@constraint(model, beta_wt[featurelookup["wt_queuecong"]] == 0)
	#@constraint(model, beta_wt[featurelookup["wt_avglocalcong"]] == 0)
	#@constraint(model, beta_wt[featurelookup["wt_avghyperlocalcong"]] == 0)
	#@constraint(model, beta_wt[featurelookup["wt_maxhyperlocalcong"]] == 0)
	#@constraint(model, beta_pwt[featurelookup["pwt_existingcong"]] == 0)

	#Isolate training
	setbeta_wt, setbeta_mp, setbeta_pw, setbeta_pwt = setbetas
	if trainassignment_flag == 1 
		@constraint(model, force_wt[f = congestionfeatures_wt], beta_wt[f] == 0)
		@constraint(model, force_pwt[f = congestionfeatures_pwt], beta_pwt[f] == 0)
	elseif traincongestion_flag == 1
		@constraint(model, force_wt[f = setdiff(featuresforprediction_wt, congestionfeatures_wt)], beta_wt[f] == setbeta_wt[f])
		@constraint(model, force_mp[f = featuresforprediction_mp], beta_mp[f] == setbeta_mp[f])
		@constraint(model, force_pw[f = featuresforprediction_pw], beta_pw[f] == setbeta_pw[f])
		@constraint(model, force_pwt[f = setdiff(featuresforprediction_pwt, congestionfeatures_pwt)], beta_pwt[f] == setbeta_pwt[f])
	end

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

	return value.(beta_wt), value.(beta_mp), value.(beta_pw), value.(beta_pwt), value.(pred), [value(beta_m), value(beta_p), value(beta_w), value(beta_t)]

end

#-----------------------------------------------------------------------------------#

function linearregressionmodel_compare(sp_workstations, sp_pods, sp_orders, sp_times, actualobj, probleminstances, numsubproblems, featurelookup, featuresforprediction_wt, featuresforprediction_mp, featuresforprediction_pw, featuresforprediction_pwt, w_features, t_features, m_features, p_features, wt_features, mp_features, pw_features, pwt_features, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums, setbetas, traincongestion_flag, trainassignment_flag, weight)

	#consistency_flag, x, y, z, v, q, actualobj = 1, length(workstationgrouplist), x_k, y_k, z_k, v_k, q_k, obj_k

	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", 60*60) # iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", 1)
	set_optimizer_attribute(model, "IntegralityFocus", 1)

	#Variables
	@variable(model, pred[i = probleminstances, k = 1:numsubproblems[i]])
	@variable(model, err[i = probleminstances, k = 1:numsubproblems[i]])
	@variable(model, beta_wt[f = wt_features])
	@variable(model, beta_mp[f = mp_features])
	@variable(model, beta_pw[f = pw_features])
	@variable(model, beta_pwt[f = pwt_features])

	@variable(model, beta_t)
	@variable(model, beta_p)
	@variable(model, beta_m)
	@variable(model, beta_w)

	@variable(model, alpha_t[f = t_featnums])
	@variable(model, alpha_p[f = p_featnums])
	@variable(model, alpha_m[f = m_featnums])
	@variable(model, alpha_w[f = w_featnums])

	@variable(model, nonzero_wt[f = featuresforprediction_wt], Bin)
	@variable(model, nonzero_mp[f = featuresforprediction_mp], Bin)
	@variable(model, nonzero_pw[f = featuresforprediction_pw], Bin)
	@variable(model, nonzero_pwt[f = featuresforprediction_pwt], Bin)
	@variable(model, chooseone, Bin)

	#Objective
	@objective(model, Min, sum(sum(weight[i,k] * err[i,k] for k in 1:numsubproblems[i]) for i in probleminstances))

	#Prediction calculation
	@constraint(model, predictioncalc[i = probleminstances, k = 1:numsubproblems[i]], pred[i,k] == 
		#Workstation-time
		sum(sum(sum(beta_wt[f] * wt_features[i,f,w,t] for f in wt_featnums) for w in sp_workstations[i,k]) for t in sp_times[i,k])
		+ sum(alpha_w[f] * w_features[i,f,w] for f in w_featnums)  
		+ sum(alpha_t[f] * t_features[i,f,t] for f in t_featnums) 
		#Order-pod
		+ sum(sum(sum(beta_mp[f] * mp_features[i,f,m,p] for f in mp_featnums) for m in sp_orders[i,k]) for p in sp_pods[i,k]) 
		+ sum(alpha_m[f] * m_features[i,f,m] for f in m_featnums)  
		+ sum(alpha_p[f] * p_features[i,f,p] for f in p_featnums) 
		#Pod-workstation
		+ sum(sum(sum(beta_pw[f] * pw_features[i,f,p,w] for f in pw_featnums) for w in sp_workstations[i,k]) for p in sp_pods[i,k]) 
		#Pod-workstation-time
		+ sum(sum(sum(sum(beta_pwt[f] * pwt_features[i,f,p,w,t] for f in pwt_featnums) for w in sp_workstations[i,k]) for p in sp_pods[i,k]) for t in sp_times[i,k])  
		#Bases
		+ sum(beta_t for t in sp_times[i,k])  
		+ sum(beta_p for t in sp_pods[i,k])
		+ sum(beta_m for t in sp_orders[i,k])
		+ sum(beta_w for t in sp_workstations[i,k])
		)

	#Prediction error
	#@constraint(model, prederror1[i = probleminstances, k = 1:numsubproblems[i]], err[i,k] >= pred[i,k] - actualobj[i,k])
	#@constraint(model, prederror2[i = probleminstances, k = 1:numsubproblems[i]], err[i,k] >= actualobj[i,k] - pred[i,k])
	for i in probleminstances, k in 1:numsubproblems[i]
		if actualobj[i,k] > 1e-4
			@constraint(model, err[i,k] >= pred[i,k] - actualobj[i,k])
			@constraint(model, err[i,k] >= actualobj[i,k] - pred[i,k])
		else
			@constraint(model, err[i,k] == 0)
		end
	end

	#Consistency / interpretability of predictions
	#@constraint(model, otherinstancesareworse[i = probleminstances, k = 1:numsubproblems[i]], pred[i, k] <= pred[i, k_star])
	@constraint(model, nonzerofeatures_wt_1[f in featuresforprediction_wt], beta_wt[f] <= 10000 * nonzero_wt[f])
	@constraint(model, nonzerofeatures_wt_2[f in featuresforprediction_wt], beta_wt[f] >= -10000 * nonzero_wt[f])
	@constraint(model, nonzerofeatures_mp_1[f in featuresforprediction_mp], beta_mp[f] <= 10000 * nonzero_mp[f])
	@constraint(model, nonzerofeatures_mp_2[f in featuresforprediction_mp], beta_mp[f] >= -10000 * nonzero_mp[f])
	@constraint(model, nonzerofeatures_pw_1[f in featuresforprediction_pw], beta_pw[f] <= 10000 * nonzero_pw[f])
	@constraint(model, nonzerofeatures_pw_2[f in featuresforprediction_pw], beta_pw[f] >= -10000 * nonzero_pw[f])
	@constraint(model, nonzerofeatures_pwt_1[f in featuresforprediction_pwt], beta_pwt[f] <= 10000 * nonzero_pwt[f])
	@constraint(model, nonzerofeatures_pwt_2[f in featuresforprediction_pwt], beta_pwt[f] >= -10000 * nonzero_pwt[f])
	@constraint(model, nonzero_mp[featurelookup["mp_itemoverlap"]] + nonzero_mp[featurelookup["mp_itemoverlappct"]] <= 1)
	@constraint(model, nonzero_mp[featurelookup["m_ordersize"]] + nonzero_mp[featurelookup["m_oneitemflag"]] <= 1)
	@constraint(model, nonzero_mp[featurelookup["mp_overlapitemavginv"]] + nonzero_mp[featurelookup["mp_overlapitemavgaltpods"]] <= 1)
	@constraint(model, nonzero_pw[featurelookup["pw_distance"]] + nonzero_pw[featurelookup["pw_xsteps"]] <= 1)
	@constraint(model, nonzero_pw[featurelookup["pw_distance"]] + nonzero_pw[featurelookup["pw_ysteps"]] <= 1)
	@constraint(model, nonzero_wt[featurelookup["wt_avglocalcong"]] + nonzero_wt[featurelookup["wt_maxlocalcong"]] <= 1)
	@constraint(model, nonzero_wt[featurelookup["wt_avghyperlocalcong"]] + nonzero_wt[featurelookup["wt_maxhyperlocalcong"]] <= 1)
	@constraint(model, nonzero_wt[featurelookup["w_centrality"]] + nonzero_pwt[featurelookup["w_centrality"]] <= 1)

	@constraint(model, beta_mp[featurelookup["p_pctusefulitems"]] == 0)
	@constraint(model, beta_pw[featurelookup["p_pctusefulitems"]] == 0)
	@constraint(model, beta_pwt[featurelookup["p_avgdisttostations"]] == 0)
	@constraint(model, beta_pw[featurelookup["p_avgdisttostations"]] == 0)
	@constraint(model, beta_pw[featurelookup["w_centrality"]] == 0)

	@constraint(model, beta_pwt[featurelookup["t_horizonrelativetime"]] == 0)
	@constraint(model, beta_wt[featurelookup["t_horizonrelativetime"]] == 0)
	@constraint(model, beta_wt[featurelookup["t_timeofday"]] == 0)
	@constraint(model, beta_pwt[featurelookup["t_timeofday"]] == 0)

	@constraint(model, beta_m == 0)
	@constraint(model, beta_p == 0)

	#@constraint(model, beta_wt[featurelookup["wt_maxlocalcong"]] == 0)
	#@constraint(model, beta_wt[featurelookup["wt_queuecong"]] == 0)
	#@constraint(model, beta_wt[featurelookup["wt_avglocalcong"]] == 0)
	#@constraint(model, beta_wt[featurelookup["wt_avghyperlocalcong"]] == 0)
	#@constraint(model, beta_wt[featurelookup["wt_maxhyperlocalcong"]] == 0)
	#@constraint(model, beta_pwt[featurelookup["pwt_existingcong"]] == 0)

	#Isolate training
	setbeta_wt, setbeta_mp, setbeta_pw, setbeta_pwt = setbetas
	if trainassignment_flag == 1 
		@constraint(model, force_wt[f = congestionfeatures_wt], beta_wt[f] == 0)
		@constraint(model, force_pwt[f = congestionfeatures_pwt], beta_pwt[f] == 0)
	elseif traincongestion_flag == 1
		@constraint(model, force_wt[f = setdiff(featuresforprediction_wt, congestionfeatures_wt)], beta_wt[f] == setbeta_wt[f])
		@constraint(model, force_mp[f = featuresforprediction_mp], beta_mp[f] == setbeta_mp[f])
		@constraint(model, force_pw[f = featuresforprediction_pw], beta_pw[f] == setbeta_pw[f])
		@constraint(model, force_pwt[f = setdiff(featuresforprediction_pwt, congestionfeatures_pwt)], beta_pwt[f] == setbeta_pwt[f])
	end

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

	return value.(beta_wt), value.(beta_mp), value.(beta_pw), value.(beta_pwt), value.(pred), [value(beta_m), value(beta_p), value(beta_w), value(beta_t)], value.(alpha_m), value.(alpha_p), value.(alpha_w), value.(alpha_t)

end

#-----------------------------------------------------------------------------------#
