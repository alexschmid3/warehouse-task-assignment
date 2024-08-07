
using JuMP, Gurobi, Plots, Random, CSV, DataFrames, Statistics, Dates, HDF5, LinearAlgebra, FileIO, JLD2, SparseArrays 

include("scripts/training/mlmodel/trainmodel.jl")

warehouse_id = parse(Int, ARGS[1])
data_pass = 0

savemodelfilename_pw = string("models/timedisc/mlmodel_wh", warehouse_id, "_pass", data_pass+1,"_pw.jld2")
savemodelfilename_nowt = string("models/timedisc/mlmodel_wh", warehouse_id, "_pass", data_pass+1,"_nowt.jld2")
savemodelfilename_full = string("models/timedisc/mlmodel_wh", warehouse_id, "_pass", data_pass+1,"_full.jld2")
savemodelfilename_intercept = string("models/timedisc/mlmodel_wh", warehouse_id, "_pass", data_pass+1,"_intercept.jld2")
#savemodelfilename_pw = string("models/newpaper/mlmodel_wh", warehouse_id, "_pass", data_pass+1,"_pw.jld2")
#savemodelfilename_nowt = string("models/newpaper/mlmodel_wh", warehouse_id, "_pass", data_pass+1,"_nowt.jld2")
#savemodelfilename_full = string("models/newpaper/mlmodel_wh", warehouse_id, "_pass", data_pass+1,"_full.jld2")
#savemodelfilename_intercept = string("models/newpaper/mlmodel_wh", warehouse_id, "_pass", data_pass+1,"_intercept.jld2")

if !(isdir("models/timedisc/"))
	mkdir("models/timedisc/")
end

#Get relevant training data
trainingfolder = string("trainingdata/timedisc_wh", warehouse_id,"/dynamic/")
unfiltered_filelist = readdir(trainingfolder)
filelist = []
for filename in unfiltered_filelist
	if occursin(string("wh", warehouse_id, "_pass", data_pass), filename)
		push!(filelist, filename)
	end
end

#Test-train split
Random.seed!(123)
numtrain = convert(Int64, round(0.8*length(filelist),digits=0))
trainfilelist = filelist[randperm(length(filelist))[1:numtrain]]
testfilelist = setdiff(filelist, trainfilelist)

#Targets
targetnumpods = 24
targetnumorders = 20
tstep = [10,20,30,40,50,60,120][warehouse_id]

allfeatures, featurenames, featurelookup, featuresforprediction_wt, featuresforprediction_mp, featuresforprediction_pw, featuresforprediction_pwt, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums = featurebank()

#Split features
congestionfeatures_wt = [featurelookup["wt_queuecong"], featurelookup["wt_avglocalcong"], featurelookup["wt_maxlocalcong"], featurelookup["wt_maxhyperlocalcong"], featurelookup["wt_avghyperlocalcong"]]
congestionfeatures_pwt = [featurelookup["pwt_existingcong"]]

#=
#Get initial model (no congestion)
traincongestion_flag = 0
trainassignment_flag = 1
#mlmodel_nocong = traindynamicsubproblemmodel(trainfilelist, traincongestion_flag, trainassignment_flag, (0,0,0,0))
mlmodel_nocong = traindynamicsubproblemmodel(filelist, traincongestion_flag, trainassignment_flag, (0,0,0,0))

#Train congestion only
traincongestion_flag = 1
trainassignment_flag = 0
#mlmodel_cong = traindynamicsubproblemmodel(trainfilelist, traincongestion_flag, trainassignment_flag, (mlmodel_nocong.beta_wt, mlmodel_nocong.beta_mp, mlmodel_nocong.beta_pw, mlmodel_nocong.beta_pwt))
mlmodel_cong = traindynamicsubproblemmodel(filelist, traincongestion_flag, trainassignment_flag, (mlmodel_nocong.beta_wt, mlmodel_nocong.beta_mp, mlmodel_nocong.beta_pw, mlmodel_nocong.beta_pwt))

println("-------------------------------------------------------------------------------------")

println("TEST")

#Test set
#testpreds, actualobj, testinstances, numsubproblems = predictdynamicsubproblemmodel(testfilelist, mlmodel_cong)
#printmodelmetrics(testpreds, actualobj, testinstances, numsubproblems)

println("-------------------------------------------------------------------------------------")

#Save model
save(savemodelfilename_pw, Dict("beta_wt" => mlmodel_cong.beta_wt, "beta_mp" => mlmodel_cong.beta_mp, "beta_pw" => mlmodel_cong.beta_pw, "beta_pwt" => mlmodel_cong.beta_pwt, "shifts" => mlmodel_cong.shifts) )

println("Piece-wise ML model saved")
=#

println("-------------------------------------------------------------------------------------")

#Train model on all data
traincongestion_flag = 0
trainassignment_flag = 0
mlmodel_nowt, mlmodel_full, mlmodel_intercept = traindynamicsubproblemmodel(filelist, traincongestion_flag, trainassignment_flag, (0,0,0,0))

println("-------------------------------------------------------------------------------------")

println("TEST")

#Test set
#testpreds, actualobj, testinstances, numsubproblems = predictdynamicsubproblemmodel(testfilelist, mlmodel_full)
#printmodelmetrics(testpreds, actualobj, testinstances, numsubproblems)

println("-------------------------------------------------------------------------------------")

function convertdenseaxistodict(denseaxisarray)
	newdict = Dict()
	for ky in keys(denseaxisarray)
		newdict[ky[1]] = denseaxisarray[ky[1]]
	end
	return newdict
end

#Save model
save(savemodelfilename_nowt, Dict("beta_wt" => convertdenseaxistodict(mlmodel_nowt.beta_wt), "beta_mp" => convertdenseaxistodict(mlmodel_nowt.beta_mp), "beta_pw" => convertdenseaxistodict(mlmodel_nowt.beta_pw), "beta_pwt" => convertdenseaxistodict(mlmodel_nowt.beta_pwt), "shifts" => mlmodel_nowt.shifts) )
#println("savemodelfilename_nowt = ")
#println(Dict("beta_wt" => mlmodel_nowt.beta_wt, "beta_mp" => mlmodel_nowt.beta_mp, "beta_pw" => mlmodel_nowt.beta_pw, "beta_pwt" => mlmodel_nowt.beta_pwt, "shifts" => mlmodel_nowt.shifts) )
save(savemodelfilename_full, Dict("beta_wt" => convertdenseaxistodict(mlmodel_full.beta_wt), "beta_mp" => convertdenseaxistodict(mlmodel_full.beta_mp), "beta_pw" => convertdenseaxistodict(mlmodel_full.beta_pw), "beta_pwt" => convertdenseaxistodict(mlmodel_full.beta_pwt), "shifts" => mlmodel_full.shifts) )
#println("savemodelfilename_full = ")
#println(Dict("beta_wt" => convertdenseaxistodict(mlmodel_nowt.beta_wt), "beta_mp" => convertdenseaxistodict(mlmodel_nowt.beta_mp), "beta_pw" => convertdenseaxistodict(mlmodel_nowt.beta_pw), "beta_pwt" => convertdenseaxistodict(mlmodel_nowt.beta_pwt), "shifts" => mlmodel_nowt.shifts) )
save(savemodelfilename_intercept, Dict("beta_wt" => convertdenseaxistodict(mlmodel_intercept.beta_wt), "beta_mp" => convertdenseaxistodict(mlmodel_intercept.beta_mp), "beta_pw" => convertdenseaxistodict(mlmodel_intercept.beta_pw), "beta_pwt" => convertdenseaxistodict(mlmodel_intercept.beta_pwt), "shifts" => mlmodel_intercept.shifts) )

println("Full ML model saved")

