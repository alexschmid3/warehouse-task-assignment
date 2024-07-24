
using JuMP, Gurobi, Plots, Random, CSV, DataFrames, Statistics, Dates, HDF5, LinearAlgebra, FileIO, JLD2, SparseArrays 

include("scripts/training/staticmlmodel/trainstaticmodel.jl")
include("scripts/training/staticmlmodel/staticmodels.jl")

warehouse_id = 6 #parse(Int, ARGS[1])
data_pass = 0
subproblemsperinstance = 1000
seed = 456
Random.seed!(seed)

savemodelfilename_obj = string("models/staticmlmodel_wh", warehouse_id, "_pass", data_pass+1,"_sp", subproblemsperinstance, "_obj.jld2")
savemodelfilename_imp = string("models/staticmlmodel_wh", warehouse_id, "_pass", data_pass+1,"_sp", subproblemsperinstance, "_imp.jld2")

#Get relevant training data
trainingfolder = string("trainingdata/cluster/mainmodel_wh",warehouse_id,"/static/")
unfiltered_filelist = readdir(trainingfolder)
filelist = []
for filename in unfiltered_filelist
	if occursin(string("wh", warehouse_id, "_pass", data_pass), filename)
		push!(filelist, filename)
	end
end

#Test-train split
numtrain = convert(Int64, round(0.8*length(filelist),digits=0))
trainfilelist = filelist[randperm(length(filelist))[1:numtrain]]
testfilelist = setdiff(filelist, trainfilelist)

#Get features
featurenames, featurenumlookup, numfeatures, trash = featurebank(trainfilelist[1])

#Find linear regression model(s)
maxselected = 12
mlmodel_obj = trainstaticsubproblemmodel_obj(filelist, maxselected)

println("------------------------------------------")

#Test predictions
#testpredictions, actualreoptobj, actualimp, originalobj, numtestsubproblems = predictdynamicsubproblemmodel(testfilelist, mlmodel_obj)
#printmodelmetrics(testpredictions, actualreoptobj, actualimp, originalobj, numtestsubproblems)

save(savemodelfilename_obj, Dict("beta" => mlmodel_obj.weights))

println("ML model saved")

println("------------------------------------------------------------------------------")

#mlmodel_imp = trainstaticsubproblemmodel_imp(filelist, maxselected)

#println("------------------------------------------")

#Test predictions
#testpredictions, actualreoptobj, actualimp, originalobj, numtestsubproblems = predictdynamicsubproblemmodel(testfilelist, mlmodel_imp)
#printmodelmetrics(testpredictions, actualreoptobj, actualimp, originalobj, numtestsubproblems)

#println("------------------------------------------")

#Save model

#save(savemodelfilename_imp, Dict("beta" => mlmodel_imp.weights))

#println("ML model saved")
