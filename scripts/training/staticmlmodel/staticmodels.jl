
function linearregressionmodel_imp(featurevalues, numsubproblems, actualimp, maxselected)

	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", 60*60*2) # iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", 1)
	set_optimizer_attribute(model, "IntegralityFocus", 1)

	#Variables
	@variable(model, pred[sp = 1:numsubproblems])
	@variable(model, err[sp = 1:numsubproblems])
	@variable(model, beta[f = 1:numfeatures])
	@variable(model, featselect[f = 1:numfeatures+1], Bin)

	#Objective
	@objective(model, Min, sum(err[sp] for sp in 1:numsubproblems))

	#Prediction calculation
	@constraint(model, predictioncalc[sp = 1:numsubproblems], pred[sp] == sum(beta[f] * featurevalues[sp,f] for f in 1:numfeatures) )

	#Prediction error
	@constraint(model, finderror1[sp = 1:numsubproblems], err[sp] >= pred[sp] - actualimp[sp])
	@constraint(model, finderror2[sp = 1:numsubproblems], err[sp] >= actualimp[sp] - pred[sp])

	#Configure choice variables
	@constraint(model, nonzerofeatures1[f in 1:numfeatures], beta[f] <= 10000 * featselect[f])
	@constraint(model, nonzerofeatures2[f in 1:numfeatures], beta[f] >= -10000 * featselect[f])

	#Sparsity
	@constraint(model, sum(featselect[f] for f in 1:numfeatures) <= maxselected)

	#Selection constraints
	@constraint(model, featselect[featurenumlookup["avg_centrality"]] + featselect[featurenumlookup["sum_centrality"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_horizonrelativetime"]] + featselect[featurenumlookup["sum_horizonrelativetime"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_avgdisttostations"]] + featselect[featurenumlookup["sum_avgdisttostations"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_ordersize"]] + featselect[featurenumlookup["sum_ordersize"]] + featselect[featurenumlookup["pct_oneitemflag"]] + featselect[featurenumlookup["num_oneitemflag"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_itemoverlappct"]] + featselect[featurenumlookup["sum_itemoverlappct"]] + featselect[featurenumlookup["avg_itemoverlap"]] + featselect[featurenumlookup["sum_itemoverlap"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_overlapitemavginv"]] + featselect[featurenumlookup["sum_overlapitemavginv"]] + featselect[featurenumlookup["avg_overlapitemavgaltpods"]] + featselect[featurenumlookup["sum_overlapitemavgaltpods"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_distance"]] + featselect[featurenumlookup["sum_distance"]] + featselect[featurenumlookup["avg_xsteps"]] + featselect[featurenumlookup["sum_xsteps"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_distance"]] + featselect[featurenumlookup["sum_distance"]] + featselect[featurenumlookup["avg_ysteps"]] + featselect[featurenumlookup["sum_ysteps"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_existingcong"]] + featselect[featurenumlookup["sum_existingcong"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_queuecong"]] + featselect[featurenumlookup["sum_queuecong"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_avglocalcong"]] + featselect[featurenumlookup["sum_avglocalcong"]] + featselect[featurenumlookup["avg_maxlocalcong"]] + featselect[featurenumlookup["sum_maxlocalcong"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_avghyperlocalcong"]] + featselect[featurenumlookup["sum_avghyperlocalcong"]] + featselect[featurenumlookup["avg_maxhyperlocalcong"]] + featselect[featurenumlookup["sum_maxhyperlocalcong"]] <= 1)
	@constraint(model, featselect[featurenumlookup["stationmissingthroughput"]] + featselect[featurenumlookup["stationidlepct"]] <= 1)
	
	#Don't pick
	@constraint(model, featselect[featurenumlookup["num_orders"]] == 0)
	@constraint(model, featselect[featurenumlookup["num_workstations"]] == 0)
	@constraint(model, featselect[featurenumlookup["num_pods"]] == 0)
	@constraint(model, featselect[featurenumlookup["horizon"]] == 0)
	@constraint(model, featselect[featurenumlookup["avgpicksperpod"]] == 0)
	@constraint(model, featselect[featurenumlookup["podsperitem"]] == 0)

	#====================================================#

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

	return value.(beta), value.(pred)

end

#-----------------------------------------------------------------------------------#

function linearregressionmodel_obj(featurevalues, numsubproblems, actualreoptobj, maxselected)

	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", 60*60*2) # iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", 1)
	set_optimizer_attribute(model, "IntegralityFocus", 1)

	#Variables
	@variable(model, pred[sp = 1:numsubproblems])
	@variable(model, err[sp = 1:numsubproblems])
	@variable(model, beta[f = 1:numfeatures+1])
	@variable(model, featselect[f = 1:numfeatures+1], Bin)

	#Objective
	@objective(model, Min, sum(err[sp] for sp in 1:numsubproblems))

	#Prediction calculation
	@constraint(model, predictioncalc[sp = 1:numsubproblems], pred[sp] == sum(beta[f] * featurevalues[sp,f] for f in 1:numfeatures+1) )

	#Prediction error
	@constraint(model, finderror1[sp = 1:numsubproblems], err[sp] >= pred[sp] - actualreoptobj[sp])
	@constraint(model, finderror2[sp = 1:numsubproblems], err[sp] >= actualreoptobj[sp] - pred[sp])

	#Configure choice variables
	@constraint(model, nonzerofeatures1[f in 1:numfeatures+1], beta[f] <= 10000 * featselect[f])
	@constraint(model, nonzerofeatures2[f in 1:numfeatures+1], beta[f] >= -10000 * featselect[f])

	#Sparsity
	@constraint(model, sum(featselect[f] for f in 1:numfeatures+1) <= maxselected)

	#Selection constraints
	@constraint(model, featselect[featurenumlookup["avg_centrality"]] + featselect[featurenumlookup["sum_centrality"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_horizonrelativetime"]] + featselect[featurenumlookup["sum_horizonrelativetime"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_avgdisttostations"]] + featselect[featurenumlookup["sum_avgdisttostations"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_ordersize"]] + featselect[featurenumlookup["sum_ordersize"]] + featselect[featurenumlookup["pct_oneitemflag"]] + featselect[featurenumlookup["num_oneitemflag"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_itemoverlappct"]] + featselect[featurenumlookup["sum_itemoverlappct"]] + featselect[featurenumlookup["avg_itemoverlap"]] + featselect[featurenumlookup["sum_itemoverlap"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_overlapitemavginv"]] + featselect[featurenumlookup["sum_overlapitemavginv"]] + featselect[featurenumlookup["avg_overlapitemavgaltpods"]] + featselect[featurenumlookup["sum_overlapitemavgaltpods"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_distance"]] + featselect[featurenumlookup["sum_distance"]] + featselect[featurenumlookup["avg_xsteps"]] + featselect[featurenumlookup["sum_xsteps"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_distance"]] + featselect[featurenumlookup["sum_distance"]] + featselect[featurenumlookup["avg_ysteps"]] + featselect[featurenumlookup["sum_ysteps"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_existingcong"]] + featselect[featurenumlookup["sum_existingcong"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_queuecong"]] + featselect[featurenumlookup["sum_queuecong"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_avglocalcong"]] + featselect[featurenumlookup["sum_avglocalcong"]] + featselect[featurenumlookup["avg_maxlocalcong"]] + featselect[featurenumlookup["sum_maxlocalcong"]] <= 1)
	@constraint(model, featselect[featurenumlookup["avg_avghyperlocalcong"]] + featselect[featurenumlookup["sum_avghyperlocalcong"]] + featselect[featurenumlookup["avg_maxhyperlocalcong"]] + featselect[featurenumlookup["sum_maxhyperlocalcong"]] <= 1)
	@constraint(model, featselect[featurenumlookup["stationmissingthroughput"]] + featselect[featurenumlookup["stationidlepct"]] <= 1)

	#Don't pick
	@constraint(model, featselect[featurenumlookup["num_orders"]] == 0)
	#@constraint(model, featselect[featurenumlookup["num_workstations"]] == 0)
	@constraint(model, featselect[featurenumlookup["num_pods"]] == 0)
	@constraint(model, featselect[featurenumlookup["old_obj"]] == 0)
	#@constraint(model, featselect[featurenumlookup["horizon"]] == 0)
	@constraint(model, featselect[featurenumlookup["avgpicksperpod"]] == 0)
	@constraint(model, featselect[featurenumlookup["podsperitem"]] == 0)

	@constraint(model, featselect[featurenumlookup["tstep"]] == 0)
	#@constraint(model, featselect[featurenumlookup["stationmissingthroughput"]] == 0)
	#@constraint(model, featselect[featurenumlookup["stationidlepct"]] == 0)
	#@constraint(model, featselect[featurenumlookup["orderslottimessopen"]] == 0)
	#@constraint(model, featselect[featurenumlookup["podcompatibilitywithunassignedorders"]] == 0)

	#@constraint(model, featselect[featurenumlookup["easyoneitemorders"]] == 0)
	#@constraint(model, featselect[featurenumlookup["totalworkedorders"]] == 0)
	#@constraint(model, featselect[featurenumlookup["stationavgordersize"]] == 0)
	#@constraint(model, featselect[featurenumlookup["stationavgorderinventory"]] == 0)
	#@constraint(model, featselect[featurenumlookup["stationpodsvisited"]] == 0)
	#@constraint(model, featselect[featurenumlookup["avgpicksperpod"]] == 0)
	@constraint(model, featselect[featurenumlookup["lsnsiterationssofar"]] == 0)
	@constraint(model, featselect[featurenumlookup["recentlyoptimized_flag"]] == 0)
	@constraint(model, featselect[featurenumlookup["totalsubproblemreoptimizations"]] == 0)
	@constraint(model, featselect[featurenumlookup["mostrecentimprovement"]] == 0)

	@constraint(model, featselect[featurenumlookup["avg_centrality"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_horizonrelativetime"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_avgdisttostations"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_ordersize"]] == 0)
	@constraint(model, featselect[featurenumlookup["pct_oneitemflag"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_itemoverlap"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_itemoverlappct"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_overlapitemavginv"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_overlapitemavgaltpods"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_distance"]] == 0)

	@constraint(model, featselect[featurenumlookup["avg_xsteps"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_ysteps"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_existingcong"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_queuecong"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_avglocalcong"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_maxlocalcong"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_avghyperlocalcong"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_maxhyperlocalcong"]] == 0)

	@constraint(model, featselect[featurenumlookup["sum_horizonrelativetime"]] == 0)

	#====================================================#

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

	return value.(beta), value.(pred)

end


#-----------------------------------------------------------------------------------#

function linearregressionmodel_compare(featurevalues, numsubproblems, actualreoptobj, maxselected)

	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", 60*60*2) # iptimelimit)
	set_optimizer_attribute(model, "OutputFlag", 1)
	set_optimizer_attribute(model, "IntegralityFocus", 1)

	#Variables
	@variable(model, pred[sp = 1:numsubproblems])
	@variable(model, err[sp = 1:numsubproblems] >= 0)
	@variable(model, beta[f = 1:numfeatures+1])
	@variable(model, featselect[f = 1:numfeatures+1], Bin)

	#Objective
	@objective(model, Min, sum(err[sp] for sp in 1:numsubproblems))

	#Prediction calculation
	for sp = 1:numsubproblems
		if actualreoptobj[sp] > 1e-4
			@constraint(model, pred[sp] == sum(beta[f] * featurevalues[sp,f] for f in 1:numfeatures+1) )
		end
	end

	#Prediction error
	@constraint(model, finderror1[sp = 1:numsubproblems], err[sp] >= pred[sp] - actualreoptobj[sp])
	@constraint(model, finderror2[sp = 1:numsubproblems], err[sp] >= actualreoptobj[sp] - pred[sp])

	#Configure choice variables
	@constraint(model, nonzerofeatures1[f in 1:numfeatures+1], beta[f] <= 10000 * featselect[f])
	@constraint(model, nonzerofeatures2[f in 1:numfeatures+1], beta[f] >= -10000 * featselect[f])

	#Selection constraints
	@constraint(model, featselect[featurenumlookup["sum_ordersize"]] + featselect[featurenumlookup["num_oneitemflag"]] <= 1)
	@constraint(model, featselect[featurenumlookup["sum_itemoverlappct"]] + featselect[featurenumlookup["sum_itemoverlap"]] <= 1)
	@constraint(model, featselect[featurenumlookup["sum_overlapitemavginv"]]+ featselect[featurenumlookup["sum_overlapitemavgaltpods"]] <= 1)
	@constraint(model, featselect[featurenumlookup["sum_distance"]] + featselect[featurenumlookup["sum_xsteps"]] <= 1)
	@constraint(model, featselect[featurenumlookup["sum_distance"]] + featselect[featurenumlookup["sum_ysteps"]] <= 1)
	@constraint(model, featselect[featurenumlookup["sum_avglocalcong"]] + featselect[featurenumlookup["sum_maxlocalcong"]] <= 1)
	@constraint(model, featselect[featurenumlookup["sum_avghyperlocalcong"]] + featselect[featurenumlookup["sum_maxhyperlocalcong"]] <= 1)

	#Don't pick
	@constraint(model, featselect[featurenumlookup["num_orders"]] == 0)
	#@constraint(model, featselect[featurenumlookup["num_workstations"]] == 0)
	@constraint(model, featselect[featurenumlookup["num_pods"]] == 0)
	@constraint(model, featselect[featurenumlookup["old_obj"]] == 0)
	#@constraint(model, featselect[featurenumlookup["horizon"]] == 0)
	@constraint(model, featselect[featurenumlookup["avgpicksperpod"]] == 0)
	@constraint(model, featselect[featurenumlookup["podsperitem"]] == 0)

	@constraint(model, featselect[featurenumlookup["tstep"]] == 0)
	@constraint(model, featselect[featurenumlookup["stationmissingthroughput"]] == 0)
	@constraint(model, featselect[featurenumlookup["stationidlepct"]] == 0)
	@constraint(model, featselect[featurenumlookup["orderslottimessopen"]] == 0)
	@constraint(model, featselect[featurenumlookup["podcompatibilitywithunassignedorders"]] == 0)

	@constraint(model, featselect[featurenumlookup["easyoneitemorders"]] == 0)
	@constraint(model, featselect[featurenumlookup["totalworkedorders"]] == 0)
	@constraint(model, featselect[featurenumlookup["stationavgordersize"]] == 0)
	@constraint(model, featselect[featurenumlookup["stationavgorderinventory"]] == 0)
	@constraint(model, featselect[featurenumlookup["stationpodsvisited"]] == 0)
	@constraint(model, featselect[featurenumlookup["avgpicksperpod"]] == 0)
	@constraint(model, featselect[featurenumlookup["lsnsiterationssofar"]] == 0)
	@constraint(model, featselect[featurenumlookup["recentlyoptimized_flag"]] == 0)
	@constraint(model, featselect[featurenumlookup["totalsubproblemreoptimizations"]] == 0)
	@constraint(model, featselect[featurenumlookup["mostrecentimprovement"]] == 0)

	@constraint(model, featselect[featurenumlookup["avg_centrality"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_horizonrelativetime"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_avgdisttostations"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_ordersize"]] == 0)
	@constraint(model, featselect[featurenumlookup["pct_oneitemflag"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_itemoverlap"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_itemoverlappct"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_overlapitemavginv"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_overlapitemavgaltpods"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_distance"]] == 0)

	@constraint(model, featselect[featurenumlookup["avg_xsteps"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_ysteps"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_existingcong"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_queuecong"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_avglocalcong"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_maxlocalcong"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_avghyperlocalcong"]] == 0)
	@constraint(model, featselect[featurenumlookup["avg_maxhyperlocalcong"]] == 0)

	@constraint(model, featselect[featurenumlookup["sum_horizonrelativetime"]] == 0)

	#====================================================#

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

	return value.(beta), value.(pred)

end