
function featurebank(firstfile)

	data = CSV.read(string(trainingfolder, firstfile), DataFrame)
	headerrow = names(data)
	featurenames = [headerrow[f] for f in 8:length(headerrow)]
	featurenumlookup = Dict()
	for f in 1:length(featurenames)
		featurenumlookup[featurenames[f]] = f
	end
	numfeatures = length(featurenames) - 5 #ignore the last five columns, which are labels, not features
	subproblemsperinstance = size(data)[1]
	df_width = size(data)[2]

	return featurenames, featurenumlookup, numfeatures, subproblemsperinstance

end

#-----------------------------------------------------------------------------------#

function readtrainingfiles(trainfilelist)
	
	featurevalues = Array{Float64}(undef, 0, numfeatures+1)
	actualreoptobj, actualimp = Array{Float64}(undef, 0, 1), Array{Float64}(undef, 0, 1)
	for file in trainfilelist 

		data = CSV.read(string(trainingfolder, file), DataFrame)
		
		run_id, instance_id = convert(Int, data[1,1]), convert(Int, data[1,2])
		dynamicfile = string(dynamictrainingfolder,"features_wh", warehouse_id,"_pass", data_pass,"_instance", instance_id,"_run", run_id,".jld2")
        println(dynamicfile)
        dynamicdata = load(dynamicfile)
		
		#Get features of all subproblems
		for sp in 1:100 #size(data)[1]
			#if data[sp,67] >= -1e-4
				featurevalues = vcat(featurevalues, Transpose(collect(data[sp,8:7+numfeatures+1])))
				actualreoptobj = vcat(actualreoptobj, data[sp,66])
				actualimp = vcat(actualimp, data[sp,67])
				spdata = dynamicdata["subproblemdata"][sp]
				spindex = size(featurevalues)[1]
				featurevalues[spindex,featurenumlookup["horizon"]] = (spdata["sp_tend"] - spdata["sp_tstart"]) / 30 + 1
				featurevalues[spindex,featurenumlookup["num_workstations"]] = length(spdata["sp_workstations"])
			#end
		end
	end
	numsubproblems = size(featurevalues)[1]
	replace!(featurevalues, NaN=>0)
	replace!(featurevalues, Inf=>0)
	replace!(featurevalues, -Inf=>0)

	return featurevalues, numsubproblems, actualreoptobj, actualimp

end	

#-----------------------------------------------------------------------------------#

function getmodelstats(beta, preds, actuals, numfeatures, numsubproblems)

	println("------------------------------------------")
	for f in 1:numfeatures
		println(featurenames[f], " = ", value(beta[f]))
	end
	println("------------------------------------------")

	totalAE, totalSE, tss, validsps = 0, 0, 0, 0
	meanobj = mean(values(actuals))
	for sp in 1:numsubproblems
		totalAE += abs(preds[sp] - actuals[sp])
		totalSE += (preds[sp] - actuals[sp])^2
		tss += (meanobj - actuals[sp])^2
		validsps += 1
	end

	println("Training mean absolute error = ", totalAE / validsps)
	println("Training mean squared error = ", totalSE / validsps)
	println("Training R_squared = ", 1 - totalSE/tss)

end

#-----------------------------------------------------------------------------------#

function trainstaticsubproblemmodel_obj(trainfilelist, maxselected)

	featurevalues, numsubproblems, actualreoptobj, actualimp = readtrainingfiles(trainfilelist)

	include("scripts/training/staticmlmodel/staticmodels.jl")
	beta, preds = linearregressionmodel_obj(featurevalues, numsubproblems, actualreoptobj, maxselected)

	getmodelstats(beta, preds, actualreoptobj, numfeatures+1, numsubproblems)

	mlmodel = (features = featurevalues, weights=beta, trainingpredictions = preds) 

	return mlmodel

end

#-----------------------------------------------------------------------------------#

function trainstaticsubproblemmodel_imp(trainfilelist, maxselected)

	featurevalues, numsubproblems, actualreoptobj, actualimp = readtrainingfiles(trainfilelist)

	include("scripts/training/staticmlmodel/staticmodels.jl")
	beta, preds = linearregressionmodel_imp(featurevalues, numsubproblems, actualimp, maxselected)

	getmodelstats(beta, preds, actualimp, numfeatures, numsubproblems)

	mlmodel = (features = featurevalues, weights=beta, trainingpredictions = preds) 

	return mlmodel
		
end

#-----------------------------------------------------------------------------------#

function trainstaticsubproblemmodel_compare(trainfilelist, maxselected)

	featurevalues, numsubproblems, actualreoptobj, actualimp = readtrainingfiles(trainfilelist)

	include("scripts/training/staticmlmodel/staticmodels.jl")
	beta, preds = linearregressionmodel_compare(featurevalues, numsubproblems, actualreoptobj, maxselected)

	getmodelstats(beta, preds, actualreoptobj, numfeatures+1, numsubproblems)

	mlmodel = (features = featurevalues, weights=beta, trainingpredictions = preds) 

	return mlmodel

end

#-----------------------------------------------------------------------------------#
#numsubproblems, beta, numfeatures, featurevalues = numsubproblems, mlmodel.weights, numfeatures, featurevalues

#=for f in 1:numfeatures
	println(featurenames[f], " --> ", featurevalues[i,f], " * ", beta[f], " = ", featurevalues[i,f] * beta[f])
end  
=#
function calcpredictionfortestset(numsubproblems, beta, numfeatures, featurevalues)

	predictions = []

	for i in 1:numsubproblems
		pred = sum(featurevalues[i,f] * beta[f] for f in 1:numfeatures)
		push!(predictions, pred)
	end

	return predictions

end

#-----------------------------------------------------------------------------------#

function predictstaticsubproblemmodel(testfilelist, mlmodel)

	featurevalues = Array{Float64}(undef, 0, numfeatures+1)
	actualreoptobj, actualimp, originalobj = Array{Float64}(undef, 0, 1), Array{Float64}(undef, 0, 1), Array{Float64}(undef, 0, 1)
	for file in testfilelist 

		data = CSV.read(string(trainingfolder, file), DataFrame)
		
		run_id, instance_id = convert(Int, data[1,1]), convert(Int, data[1,2])

		dynamicfile = string(dynamictrainingfolder,"features_wh", warehouse_id,"_pass", data_pass,"_instance", instance_id,"_run", run_id,".jld2")
        println(dynamicfile)
        dynamicdata = load(dynamicfile)

		#Get features of all subproblems
		for sp in 1:100 #size(data)[1]
			#if data[sp,67] >= -1e-4
				featurevalues = vcat(featurevalues, Transpose(collect(data[sp,8:7+numfeatures+1])))
				actualreoptobj = vcat(actualreoptobj, data[sp,66])
				actualimp = vcat(actualimp, data[sp,67])
				originalobj = vcat(originalobj, data[sp,65])
				spdata = dynamicdata["subproblemdata"][sp]
				spindex = size(featurevalues)[1]
				featurevalues[spindex,featurenumlookup["horizon"]] = (spdata["sp_tend"] - spdata["sp_tstart"])  / 30 + 1
				featurevalues[spindex,featurenumlookup["num_workstations"]] = length(spdata["sp_workstations"])
			#end
		end
	end
	numsubproblems = size(featurevalues)[1]
	replace!(featurevalues, NaN=>0)
	replace!(featurevalues, Inf=>0)
	replace!(featurevalues, -Inf=>0)

	testpredictions = calcpredictionfortestset(numsubproblems, mlmodel.weights, numfeatures, featurevalues)

	return testpredictions, actualreoptobj, actualimp, originalobj, numsubproblems

end

#-----------------------------------------------------------------------------------#

function printmodelmetrics_static(preds, actual, actualimp, originalobj, numsubproblems)
	
	totalAE, totalSE, tss, validsps = 0, 0, 0, 0
	meanobj = mean(values(actual))
	for sp in 1:numsubproblems
		if actual[sp] >= 1e-4
			totalAE += abs(preds[sp] - actual[sp])
			totalSE += (preds[sp] - actual[sp])^2
			tss += (meanobj - actual[sp])^2
			validsps += 1
		end
	end

	println("Testing mean absolute error = ", totalAE / validsps)
	println("Testing mean squared error = ", totalSE / validsps)
	println("Testing R_squared = ", 1 - totalSE/tss)

	top10_act = [sp for sp in 1:numsubproblems if actualimp[sp] >= sort(actualimp, dims=1)[numsubproblems - convert(Int, ceil(numsubproblems*0.1))]]
	top5_act = [sp for sp in 1:numsubproblems if actualimp[sp] >= sort(actualimp, dims=1)[numsubproblems - convert(Int, ceil(numsubproblems*0.05))]]
	top1_act = [sp for sp in 1:numsubproblems if actualimp[sp] >= sort(actualimp, dims=1)[numsubproblems - convert(Int, ceil(numsubproblems*0.01))]]
	top01_act = [sp for sp in 1:numsubproblems if actualimp[sp] >= sort(actualimp, dims=1)[numsubproblems - convert(Int, ceil(numsubproblems*0.001))]]

	predimp = [preds[sp] - originalobj[sp] for sp in 1:numsubproblems]
	top10_pred = [sp for sp in 1:numsubproblems if preds[sp] - originalobj[sp] >= sort(predimp)[numsubproblems - convert(Int, ceil(numsubproblems*0.1))]]
	top5_pred = [sp for sp in 1:numsubproblems if preds[sp] - originalobj[sp] >= sort(predimp)[numsubproblems - convert(Int, ceil(numsubproblems*0.05))]]
	top1_pred = [sp for sp in 1:numsubproblems if preds[sp] - originalobj[sp] >= sort(predimp)[numsubproblems - convert(Int, ceil(numsubproblems*0.01))]]
	top01_pred = [sp for sp in 1:numsubproblems if preds[sp] - originalobj[sp] >= sort(predimp)[numsubproblems - convert(Int, ceil(numsubproblems*0.001))]]

	println("Top 10 pct identification accuracy = ", length(intersect(top10_act, top10_pred)) / length(top10_pred) )
	println("Top 5 pct identification accuracy = ", length(intersect(top5_act, top5_pred)) / length(top5_pred) )
	println("Top 1 pct identification accuracy = ", length(intersect(top1_act, top1_pred)) / length(top1_pred) )
	println("Top 0.1 pct identification accuracy = ", length(intersect(top01_act, top01_pred)) / length(top01_pred) )

end




