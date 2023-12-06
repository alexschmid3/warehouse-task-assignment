
function enumeratesubproblemcandidates(subproblembudget, currpartition, windows, windowidlookup, currsol)

    spselectionchoices = []
	for sp_index in 1:subproblembudget
        sp_winid, sp_orders, sp_window, sp_pods, sp_itemson, sp_items = selectrandomsubproblem(currpartition, windows, windowidlookup, currsol, targetnumorders, targetnumpods)
        sp_workstations, sp_times = sp_window

        sp = constructsubproblem(currpartition, sp_orders, sp_window, sp_pods, sp_itemson, sp_items, currsol)

        oldobj = sum(sum(length(currsol.itempodpicklist[w,t]) for t in sp_times) for w in sp_workstations)
        featurevals = calcfeatures()
        
        spfull = (oldobj=oldobj, features=featurevals, sp=sp)
        
        push!(spselectionchoices, spfull)
    end

    return spselectionchoices

end

#-----------------------------------------------------------------------------------------------------#

function pickbestsubproblem(spselectionchoices)

    predicted_objs, old_objs = [], []
	numfeatures = 57

	for sp in spselectionchoices
		featurevalues = sp.features[5:4+numfeatures]
		pred = sum(featurevalues[f] * beta[f] for f in 1:numfeatures)
		push!(predicted_objs, pred)
		push!(old_objs, sp.old_obj)
	end

	predicted_imp = [predicted_objs[sp] - old_objs[sp] for sp in 1:length(spselectionchoices)]

	sp_index = argmax(predicted_imp)
	chosensp = spselectionchoices[sp_index]

	return chosensp.winid, chosensp.orders, chosensp.window, chosensp.window, chosensp.itemson, chosensp.items

end

#-----------------------------------------------------------------------------------------------------#

function selectlearningbenchmarksubproblem(currpartition, currsol, windows, windowscontaining, windowidlookup, windowsduring, windowsynergy, targetnumorders, targetnumpods, tabulist, lastoptimizeddifference, sp_iter)

    spselectionchoices = enumeratesubproblemcandidates(subproblembudget, currpartition, windows, windowidlookup, currsol)
    sp_winid, sp_orders, sp_window, sp_pods, sp_itemson, sp_items = pickbestsubproblem(spselectionchoices)

	return sp_winid, sp_orders, sp_window, sp_pods, sp_itemson, sp_items, tabulist, objective

end