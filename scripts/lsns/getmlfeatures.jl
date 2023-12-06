
function featurebank()

	w_features, t_features, m_features, p_features = Dict(), Dict(), Dict(), Dict()
	wt_features, mp_features, pw_features, pwt_features = Dict(), Dict(), Dict(), Dict(), Dict()
	actualobj, actualimp = Dict(), Dict()
	oldobj = Dict()
	sp_workstations, sp_pods, sp_times, sp_orders = Dict(), Dict(), Dict(), Dict(), Dict()

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

#-----------------------------------------------------------------------------#

function getinstancefeatures(currpartition, currsol)

	currentobjective = sum(sum(length(currsol.itempodpicklist[w,t]) for t in times) for w in currpartition.workstations)
	timeused = podprocesstime * sum(sum(length(currsol.podsworkedat[w,t]) for t in times) for w in currpartition.workstations) + itemprocesstime * sum(sum(length(currsol.itempodpicklist[w,t]) for t in times) for w in currpartition.workstations)
 
	instance_features = Dict("lsnsiterationsbefore" => lsnsiterationsbefore,
							"dynamicmlpass" => dynamicmlpass, 
							"greedymethod" => orderprioritization, 
							"throughpututilization" => currentobjective / (horizon * numworkstations / (podprocesstime + itemprocesstime)), 
							"timeutilization" => timeused / (horizon * numworkstations), 
							"podsperitem" => numpods * num_items_per_pod / num_unique_items , 
							"orders" => orders, 
							"itemson" => itemson, 
							"podswith" => podswith, 
							"pods" => pods, 
							"workstations" => workstations, 
							"tstep" => tstep, 
							"horizon" => horizon)

	return instance_features

end

#-----------------------------------------------------------------------------#

function getentityfeatures()
	 
	w_features = Dict("centrality" => Dict())
	t_features = Dict("horizonrelativetime" => Dict(), "timeofday" => Dict())
	p_features = Dict("avgdisttostations" => Dict(), "pctusefulitems" => Dict())
	m_features = Dict("ordersize" => Dict(), "oneitemflag" => Dict())
	i_features = Dict("iteminventory" => Dict(), "itemnumpods" => Dict())

	for w in workstations 
		w_features["centrality"][w] = sum(warehousedistance[l,w] for l in storagelocs) / length(storagelocs)
	end

	for t in extendedtimes
		t_features["horizonrelativetime"][t] = t / horizon
		t_features["timeofday"][t] = 0
	end

	for p in pods
		p_features["avgdisttostations"][p] = mean([manhattandist(podstorageloc[p],w) for w in workstations])
		p_features["pctusefulitems"][p] = length(podstartinventory[p]) / num_items_per_pod
	end

	for m in orders
		m_features["ordersize"][m] = length(itemson[m])
		m_features["oneitemflag"][m] = 2 - min(2, length(itemson[m]))
	end

	for i in items
		i_features["iteminventory"][i] = sum(inventory[i,p] for p in podswith[i])
		i_features["itemnumpods"][i] = length(podswith[i])
	end

	return w_features, t_features, p_features, m_features, i_features

end

#-----------------------------------------------------------------------------#

function getsynergyfeatures()
	
	mp_features = Dict("itemoverlap" => spzeros(length(orders), length(pods)), "itemoverlappct" => spzeros(length(orders), length(pods)), "overlapitemavginv" => spzeros(length(orders), length(pods)), "overlapitemavgaltpods" => spzeros(length(orders), length(pods)))
	pw_features = Dict("distance" => zeros(length(pods), length(workstations)), "xsteps" => zeros(length(pods), length(workstations)), "ysteps" => zeros(length(pods), length(workstations)))
	pwt_features = Dict("existingcong" => zeros(length(pods), length(workstations), length(times)))
	wt_features = Dict("queuecong" => spzeros(length(workstations), length(times)), "avglocalcong" => spzeros(length(workstations), length(times)), "maxlocalcong" => spzeros(length(workstations), length(times)), "avghyperlocalcong" => spzeros(length(workstations), length(times)), "maxhyperlocalcong" => spzeros(length(workstations), length(times)))

	#Order-pod inventory features
	itemcounter = spzeros(length(orders), length(pods))
	for m in orders, i in itemson[m], p in podswith[i]
		mp_features["itemoverlap"][m,p] += inventory[i,p]
		mp_features["overlapitemavginv"][m,p] += sum(inventory[i,p2] for p2 in podswith[i]) 
		mp_features["overlapitemavgaltpods"][m,p] += length(podswith[i])-1
		itemcounter[m,p] += 1
	end
	for (m,p) in zip(findnz(mp_features["itemoverlap"])...)
		mp_features["itemoverlappct"][m,p] = mp_features["itemoverlap"][m,p] / length(itemson[m])
		mp_features["overlapitemavginv"][m,p] = mp_features["overlapitemavginv"][m,p] / itemcounter[m,p]
		mp_features["overlapitemavgaltpods"][m,p] = mp_features["overlapitemavgaltpods"][m,p] / itemcounter[m,p]
	end

	#Pod-workstation features 
	for p in pods, w in workstations
		pw_features["distance"][p,w-last(storagelocs)] = manhattandist(podstorageloc[p],w)
		pw_features["xsteps"][p,w-last(storagelocs)] = abs(loccoords[podstorageloc[p],1] - loccoords[w,1]) 
		pw_features["ysteps"][p,w-last(storagelocs)] = abs(loccoords[podstorageloc[p],2] - loccoords[w,2])
	end

	#Pod-workstation-time congestion
	
	#for s in storagelocs[1:2], w in workstations, t in 0:tstep:horizon-arclength[s,w]
	#	a = arcs[nodes[s,t], nodes[w,t+arclength[s,w]]]
	#	relevantcongestion = congestionsignature[a] .* sum(currcong[p] for p in pods) ./ maxcong
	#	if sum(relevantcongestion) > 1e-4
	#		for p in podsstoredat[s]
	#			pwt_features["existingcong"][p,w - last(storagelocs), t/tstep + 1] = sum(relevantcongestion) / length(findnz(relevantcongestion)[1])
	#		end
	#	end
	#end

	#Congestion pre-processing
	allcong = sum(currcong[p] for p in pods)
	maxcong = zeros(length(intersections), length(0:congestiontstep:horizon))
	for i in intersections, t in 1:length(0:congestiontstep:horizon)
		maxcong[maps.mapintersectiontorow[i], t] = intersectionmaxpods[i]
	end
	normalizedcong = allcong ./ maxcong

	#Workstation-time congestion features
	kdtree = KDTree(transpose(intcoords_nn))
	for w in workstations, t in 0:tstep:horizon
		#Queue congestion
		q = maps.mapintersectiontorow[maploctointersection[w]]
		wt_features["queuecong"][w-last(storagelocs),convert(Int,t/tstep+1)] = mean([normalizedcong[q,maps.maptimetocolumn[t2]] for t2 in max(0,t-tstep):congestiontstep:t])
	
		#Local congestion
		closeints, dists = knn(kdtree, [intcoords[q][1], intcoords[q][2]], numlocalintersections, true)
		congestionnearworkstation = [sum(normalizedcong[maps.mapintersectiontorow[i+minimum(intersections)-1],maps.maptimetocolumn[t2]] for t2 in max(0,t-tstep):congestiontstep:t) for i in closeints]
		wt_features["avglocalcong"][w-last(storagelocs),convert(Int,t/tstep+1)] = mean(congestionnearworkstation)
		wt_features["maxlocalcong"][w-last(storagelocs),convert(Int,t/tstep+1)] = maximum(congestionnearworkstation)

		#Hyper local congestion
		closestints, dists = knn(kdtree, [intcoords[q][1], intcoords[q][2]], numhyperlocalintersections, true)
		congestionnearestworkstation = [sum(normalizedcong[maps.mapintersectiontorow[i+minimum(intersections)-1],maps.maptimetocolumn[t2]] for t2 in max(0,t-tstep):congestiontstep:t) for i in closestints]
		wt_features["avghyperlocalcong"][w-last(storagelocs),convert(Int,t/tstep+1)] = mean(congestionnearestworkstation)
		wt_features["maxhyperlocalcong"][w-last(storagelocs),convert(Int,t/tstep+1)] = maximum(congestionnearestworkstation)
	end

	return mp_features, pw_features, pwt_features, wt_features

end

#-----------------------------------------------------------------------------------------------------#

function getmlfeatures(mlmodelfilename)

	#ML coefficients
	#mldata = load(mlmodelfilename)
	#beta_wt, beta_mp, beta_pw, beta_pwt, beta_t, beta_w = mldata["beta_wt"], mldata["beta_mp"], mldata["beta_pw"], mldata["beta_pwt"], mldata["shifts"][4], mldata["shifts"][3]
	#beta = (wt=beta_wt, mp=beta_mp, pw=beta_pw, pwt=beta_pwt, t=beta_t, w=beta_w)
	if mlmodelfilename == "models/mlmodel_wh2.jld2"
		beta = (wt = Dict(13 => 2.5600584713264474, 2 => 0.0, 10 => -0.40109211853918064, 11 => 0.0, 12 => -0.05043682643320213, 14 => 0.0, 3 => 0.0, 1 => 0.607447336638065), mp=Dict(5 =>  0.0, 4 => -0.002176404137823992, 16 => 0.0, 7 => 0.19103193629980503, 15 => 0.08261184720499964, 6 =>  0.0, 18 => -0.005315607731556266, 17 => 0.0), t = 2.8759661851642013, w = -11.890577826920186, pw = Dict(5 => 0.0, 4 => 0.0, 21 => 0.018609090014792447, 20 => -0.012013416644662994, 1 => 0.0, 19 => 0.0), pwt = Dict(5 => -0.704737897235883, 4 => 0.0, 22 => 0.040319450361324456, 2 => 0.0, 3 => 0.0, 1 => 0.0))
		#beta = (wt = Dict(13 => 0.0, 2 => 0.0, 10 => -1.4961831312534941, 11 => 3.4490407010271285, 12 => 0.0, 14 => 1.4237052201931832, 3 => 0.0, 1 => -0.06685921041384356), mp=Dict(5 => 0.0, 4 => -0.004418310859418085, 16 => 0.0, 7 => 0.07945124286152572, 15 => 0.12974755697656148, 6 => 0.0, 18 => 0.0, 17 => -0.008203402756941756), t = 3.065438071374471, w = -0.6298225932825006, pw = Dict(5 => 0.0, 4 => 0.0, 21 => 0.024773565764949153, 20 => 6.456436723863672e-5, 1 => 0.0, 19 => 0.0), pwt = Dict(5 => -0.7131759433379721, 4 => 0.0, 22 => 0.0, 2 => 0.0, 3 => 0.0, 1 => 0.0))
	elseif mlmodelfilename == "models/mlmodel_wh1.jld2"
		beta = (wt = Dict(13 => 2.48201365992898, 2 => 0.0, 10 => -1.1974822922837314, 11 => 0.0, 12 => 1.1539562003741552, 14 => 0.0, 3 => 0.0, 1 => 0.0), mp=Dict(5 =>  0.0, 4 =>  -0.002382875732504249, 16 => 0.0, 7 => 0.0, 15 =>  0.4638863268345246, 6 => -0.14728761083001365, 18 => 0.0, 17 => -0.024265008594418363), t = 0.5208745061443563, w = -2.222707517366639, pw = Dict(5 => 0.0, 4 => 0.0, 21 =>  0.004746001059430387, 20 => -0.005054752125122445, 1 => 0.0, 19 => 0.0), pwt = Dict(5 =>  -0.07553631399052436, 4 => 0.0, 22 => 0.0, 2 => 0.0, 3 => 0.0, 1 => 0.10892834905935578))
	elseif mlmodelfilename == "models/mlmodel_wh5.jld2"
		beta = (wt = Dict(13 => 5.512869156834523, 2 => 0.0, 10 => -1.42108323426685, 11 => 0.0, 12 => 0.024233372196918224, 14 => 0.0, 3 => 0.0, 1 => 0.04495363606281663), mp=Dict(5 => 0.0, 4 => -0.016540556945368504, 16 => 0.0, 7 => 0.0, 15 => 0.22602321484690469, 6 => 0.07910873906447316, 18 => 0.0, 17 => -0.02173859857859117), t = 1.3733099495288186, w = -4.929254481374822, pw = Dict(5 => 0.0, 4 => 0.0, 21 => 0.01254264546096824, 20 => -0.005553885900396432, 1 => 0.0, 19 => 0.0), pwt = Dict(5 => 0.05092503118533457, 4 => 0.0, 22 => 0.0, 2 => 0.0, 3 => 0.0, 1 => 0.0))
	elseif mlmodelfilename == "models/mlmodel_wh6.jld2" #CHANGE
		beta = (wt = Dict(13 => 0.0, 2 => 0.0, 10 => -1.4843872816077883, 11 => 3.9807897088336603, 12 => 0.0, 14 => 0.14399618997421393, 3 => 0.0, 1 => 0.0), mp=Dict(5 => 0.0, 4 => 0.00024822623482111567, 16 => 0.0, 7 => 0.0, 15 => 0.9365828622354935, 6 => -0.06830408963527207, 18 => 0.0, 17 => -0.06361387003541773), t = 0.30144160647503204, w = -4.786219502341641, pw = Dict(5 => 0.0, 4 => 0.0, 21 => -0.0014407386621202221, 20 => -0.004702039790482948, 1 => 0.0, 19 => 0.0), pwt = Dict(5 => -0.10679737300365158, 4 => 0.0, 22 => 0.0, 2 => 0.0, 3 => 0.0, 1 => 0.20767820909893453))
	elseif mlmodelfilename == "models/mlmodel_wh3.jld2"
		beta = (wt = Dict(13 => 2.5600584713264474, 2 => 0.0, 10 => -0.40109211853918064, 11 => 0.0, 12 => -0.05043682643320213, 14 => 0.0, 3 => 0.0, 1 => 0.607447336638065), mp=Dict(5 =>  0.0, 4 => -0.002176404137823992, 16 => 0.0, 7 => 0.19103193629980503, 15 => 0.08261184720499964, 6 =>  0.0, 18 => -0.005315607731556266, 17 => 0.0), t = 2.8759661851642013, w = -11.890577826920186, pw = Dict(5 => 0.0, 4 => 0.0, 21 => 0.018609090014792447, 20 => -0.012013416644662994, 1 => 0.0, 19 => 0.0), pwt = Dict(5 => -0.704737897235883, 4 => 0.0, 22 => 0.040319450361324456, 2 => 0.0, 3 => 0.0, 1 => 0.0))
	elseif mlmodelfilename == "models/mlmodel_wh4.jld2"
		beta = (wt = Dict(13 => 5.512869156834523, 2 => 0.0, 10 => -1.42108323426685, 11 => 0.0, 12 => 0.024233372196918224, 14 => 0.0, 3 => 0.0, 1 => 0.04495363606281663), mp=Dict(5 => 0.0, 4 => -0.016540556945368504, 16 => 0.0, 7 => 0.0, 15 => 0.22602321484690469, 6 => 0.07910873906447316, 18 => 0.0, 17 => -0.02173859857859117), t = 1.3733099495288186, w = -4.929254481374822, pw = Dict(5 => 0.0, 4 => 0.0, 21 => 0.01254264546096824, 20 => -0.005553885900396432, 1 => 0.0, 19 => 0.0), pwt = Dict(5 => 0.05092503118533457, 4 => 0.0, 22 => 0.0, 2 => 0.0, 3 => 0.0, 1 => 0.0))
	end
	
	#Entity features
	allfeatures, featurenames, featurelookup, featuresforprediction_wt, featuresforprediction_mp, featuresforprediction_pw, featuresforprediction_pwt, w_featnums, t_featnums, p_featnums, m_featnums, wt_featnums, mp_featnums, pw_featnums, pwt_featnums = featurebank()
	#instance_features = getinstancefeatures(currpartition, currsol)
	w_features, t_features, p_features, m_features, i_features = getentityfeatures()
	mp_features, pw_features, pwt_features, wt_features = getsynergyfeatures()

	featureinfo = (listall=allfeatures, names=featurenames, lookup=featurelookup)
	featurenums = (w=w_featnums, t=t_featnums, p=p_featnums, m=m_featnums, wt=wt_featnums, mp=mp_featnums, pw=pw_featnums, pwt=pwt_featnums)
	featuresfor = (wt=featuresforprediction_wt, mp=featuresforprediction_mp, pw=featuresforprediction_pw, pwt=featuresforprediction_pwt)
	features = (w=w_features, t=t_features, p=p_features, m=m_features, i=i_features, mp=mp_features, pw=pw_features, pwt=pwt_features, wt=wt_features) #, currsol=currsol_features, openorders=openorders, openpods=openpods)

	return beta, features, featuresfor, featurenums, featureinfo

end

# How to organize partition features once the partitions are set?
# Debug overall method
# Kickoff 1 partition and 10 partitions on cluster by Wednesday
