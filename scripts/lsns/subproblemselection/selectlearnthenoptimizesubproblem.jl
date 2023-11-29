
include("findassignedorders.jl")

#-----------------------------------------------------------------------------------------------------#

function preprocesswindowsynergies(windows, windowidlookup, featureinfo, features, beta, featurenums)

	windowsynergy_wt, windowsynergy_w, windowsynergy_t  = Dict(), Dict(), Dict()

	for win in windows
		windowsynergy_wt[windowidlookup[win]] = sum(sum(sum(beta.wt[f] * features.wt[chop(featureinfo.names[f],head=3,tail=0)][w-last(storagelocs),convert(Int,t/tstep+1)] for f in featurenums.wt) for w in win.workstations) for t in intersect(times, win.times)) 
		windowsynergy_wt[windowidlookup[win]] += sum(sum(beta.wt[f] * features.w[chop(featureinfo.names[f],head=2,tail=0)][w] for f in featurenums.w) for w in win.workstations) 
		windowsynergy_wt[windowidlookup[win]] += sum(sum(beta.wt[f] * features.t[chop(featureinfo.names[f],head=2,tail=0)][t] for f in featurenums.t) for t in win.extendedtimes)
		windowsynergy_t[windowidlookup[win]] = sum(beta.t for t in win.times)
		windowsynergy_w[windowidlookup[win]] = sum(beta.w for w in win.workstations) 
	end

	windowsynergy = (wt=windowsynergy_wt, w=windowsynergy_w, t=windowsynergy_t)

	return windowsynergy

end

#-----------------------------------------------------------------------------------------------------#

function getassinedordersbystationtime(currpartition, currsol)

	assignedorders = Dict()
	for w in currpartition.workstations, t in times
		#Include all open orders, regardless of whether there are items picked in the window
		assignedorders[w,t] = currsol.ordersopen[w,t]

		#Include all orders that have an item picked in the window, regardless of whether they remained open for a full period of the window
		#Note this is necessary to ensure single item orders are included
		for (m,i,p) in currsol.itempodpicklist[w,t]
			assignedorders[w,t] = union(assignedorders[w,t], m)
		end
	end

	return assignedorders

end

#-----------------------------------------------------------------------------------------------------#

function constructbestsubproblem(currpartition, currsol, features, featurenums, featureinfo, windows, windowscontaining, windowidlookup, windowsduring, windowsynergy, unassignedorders, tabulist, lastoptimizeddifference, sp_iter)

	assignedorders = getassinedordersbystationtime(currpartition, currsol)
	windowids = 1:length(windows)

	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "TimeLimit", timeforsubproblemselection)
	set_optimizer_attribute(model, "OutputFlag", 0)
	set_optimizer_attribute(model, "MIPGap", 0.03)

	#Variables
	@variable(model, x[windowids], Bin)
	@variable(model, alpha[currpartition.orders], Bin)
	@variable(model, zeta[currpartition.pods], Bin)
	@variable(model, y[currpartition.orders, currpartition.pods], Bin)
	@variable(model, z[currpartition.pods, windowids], Bin)
	@variable(model, tabupenalty[1:length(tabulist)], Bin)
	@variable(model, podpenalty >= 0)
	@variable(model, itempenalty >= 0)
	@variable(model, predicted_objective)

	#Objective = maximize synergy minus current objective
	@objective(model, Max, 
		predicted_objective 
		#Current objective
		- sum(sum(sum((1-rand()/5) * length(currsol.itempodpicklist[w,t]) * x[windowidlookup[win]] for t in win.times) for w in win.workstations) for win in windows)
		#Tabu penalty 
		- 40*sum(tabupenalty[sp] for sp in 1:length(tabulist))
		- 10*podpenalty - 5*itempenalty
		- sum(lastoptimizeddifference[win] * x[win] for win in tabulist))

	@constraint(model, predicted_objective == 
		#Workstation-time synergy
		sum(windowsynergy.wt[win] * x[win] for win in windowids) 
		#Order-pod synergy
		+ sum(sum(sum(beta.mp[f] * features.mp[chop(featureinfo.names[f],head=3,tail=0)][m,p] * y[m,p] for f in featurenums.mp) for m in currpartition.orders) for p in currpartition.pods) 
		+ sum(sum(beta.mp[f] * features.m[chop(featureinfo.names[f],head=2,tail=0)][m] * alpha[m] for f in featurenums.m) for m in currpartition.orders) 
		+ sum(sum(beta.mp[f] * features.p[chop(featureinfo.names[f],head=2,tail=0)][p] * zeta[p] for f in featurenums.p) for p in currpartition.pods) 
		#Pod-workstation synergy
		+ sum(sum(sum(sum(beta.pw[f] * features.pw[chop(featureinfo.names[f],head=3,tail=0)][p,w-last(storagelocs)] * z[p,windowidlookup[win]] for f in featurenums.pw) for w in win.workstations) for win in windows) for p in currpartition.pods) 
		+ sum(sum(beta.pw[f] * features.p[chop(featureinfo.names[f],head=2,tail=0)][p] * zeta[p] for f in featurenums.p) for p in currpartition.pods) 
		+ sum(sum(sum(beta.pw[f] * features.w[chop(featureinfo.names[f],head=2,tail=0)][w] * x[windowidlookup[win]] for f in featurenums.w) for w in win.workstations) for win in windows)
		#Pod-workstation-time synergy
		+ sum(sum(sum(sum(sum(beta.pwt[f] * features.pwt[chop(featureinfo.names[f],head=4,tail=0)][p,w-last(storagelocs),convert(Int,t/tstep+1)] * z[p,windowidlookup[win]] for f in featurenums.pwt) for w in win.workstations) for p in currpartition.pods) for t in win.times) for win in windows)
		+ sum(sum(beta.pwt[f] * features.p[chop(featureinfo.names[f],head=2,tail=0)][p] * zeta[p] for f in featurenums.p) for p in currpartition.pods)
		+ sum(sum(sum(beta.pwt[f] * features.w[chop(featureinfo.names[f],head=2,tail=0)][w] * x[windowidlookup[win]] for f in featurenums.w) for w in win.workstations) for win in windows)  
		+ sum(sum(sum(beta.pwt[f] * features.t[chop(featureinfo.names[f],head=2,tail=0)][t] * x[windowidlookup[win]] for f in featurenums.t) for t in win.extendedtimes) for win in windows)
		#Time blocks
		+ sum(windowsynergy.t[win] * x[win] for win in windowids)
		+ sum(windowsynergy.w[win] * x[win] for win in windowids))

	#Force specific window to be chosen
	if windowforcingflag == 1
		@constraint(model, x[mod(sp_iter-1, length(windows))+1] == 1)
	end

	#Selection consistency
	@constraint(model, chooseonewindow, sum(x[win] for win in windowids) == 1)
	@constraint(model, consistency_mp[m in currpartition.orders, p in currpartition.pods], 2 * y[m,p] <= alpha[m] + zeta[p])
	@constraint(model, consistency_pw1[p in currpartition.pods, win in windowids], 2 * z[p,win] <= x[win] + zeta[p])
	@constraint(model, consistency_pw2[p in currpartition.pods, win in windowids], z[p,win] >= x[win] + zeta[p] - 1)

	#Order-item-pod consistency
	@constraint(model, orderitempod[m in currpartition.orders, i in itemson[m]], alpha[m] <= sum(zeta[p] for p in currpartition.podswith[i]) )

	#Order and pod inclusion
	@constraint(model, orderinclusion[w in currpartition.workstations, t in times, m in unique(assignedorders[w,t]), win in windowscontaining[w,t]], alpha[m] >= x[win] )
	@constraint(model, podinclusion[w in currpartition.workstations, t in times, p in unique(currsol.podsworkedat[w,t]), win in windowscontaining[w,t]], zeta[p] >= x[win] )
	
	#Order and pod exclusion
	@constraint(model, orderisnotthere[m in setdiff(currpartition.orders, unassignedorders)], alpha[m] <= 0) 
	for w in currpartition.workstations, t in times, m in assignedorders[w,t]
		for win in windowscontaining[w,t]
			set_normalized_coefficient(orderisnotthere[m], x[win], -1.0)
		end
	end
	@constraint(model, podisnotthere[p in currpartition.pods, (w,t) in currsol.podbusy[p], win in setdiff(windowsduring[t], windowscontaining[w,t])], z[p,win] == 0)

	#Subproblem size
	#@constraint(model, maxorders, sum(alpha[m] for m in orders) <= targetnumorders) #length(orders) / 5 )
	#@constraint(model, minorders, sum(alpha[m] for m in orders) >= minnumorders) #length(orders) / 5 )
	@constraint(model, maxpods, sum(zeta[p] for p in currpartition.pods) <= targetnumpods + podpenalty) #length(pods) / 10 )
	@constraint(model, minpods, sum(zeta[p] for p in currpartition.pods) >= minnumpods - podpenalty) #length(pods) / 10 )
	@constraint(model, maxitems, sum(length(itemson[m]) * alpha[m] for m in currpartition.orders) <= targetnumitems + itempenalty) 
	@constraint(model, minitems, sum(length(itemson[m]) * alpha[m] for m in currpartition.orders) >= minnumitems - itempenalty) 
	@constraint(model, podpenalty <= 10)
	@constraint(model, itempenalty <= 20)

	#Taboo workstations and times
	for sp in 1:length(tabulist)
		sp_window = tabulist[sp]
		@constraint(model, x[sp_window] <= tabupenalty[sp]) 
	end

	#====================================================#

	#Solve IP
	status = optimize!(model)

	if (termination_status(model) == MOI.OPTIMAL) || ((termination_status(model) == MOI.TIME_LIMIT) && (has_values(model)))
		return value(predicted_objective), x, alpha, zeta
	elseif (termination_status(model) == MOI.INFEASIBLE) #|| (termination_status(model) == MOI.INF_OR_UNBD)
		println("Infeasible SP selection!")
		return nothing
	else
		println("No solution found!")
		println(termination_status(model))
		return nothing
	end

end

#-----------------------------------------------------------------------------------------------------#

function parseoptimizationoutput(windows, currpartition, currsol, x, alpha, zeta, tabulist)

	#Get the selected entities
	sp_winid = [win for win in 1:length(windows) if value(x[win]) > 1e-4][1]
	sp_window = windows[sp_winid]
	sp_orders = [m for m in currpartition.orders if value(alpha[m]) > 1e-4]	
	sp_pods = [p for p in currpartition.pods if value(zeta[p]) > 1e-4]	

	#Find the relevant items
	assignedorders, unassignedorders_trash, numadditionalorders_trash = findassignedorders(currpartition, sp_window, currsol, 0)
	sp_itemson = Dict()
	for m in assignedorders
		sp_itemson[m] = []
		for i in itemson[m]
			if checkiteminpicklist(m, i, max(0,sp_window.tstart), min(horizon,sp_window.tend), sp_window.workstations, currsol)
				push!(sp_itemson[m], i)
			elseif (sp_window.tstart == -30) && (sum(sum(currsol.h[m,i,p,w,sp_window.tend] for p in currpartition.podswith[i]) for w in sp_window.workstations) > 0.01)
				push!(sp_itemson[m], i)
			end
		end
	end

	additionalorders = setdiff(sp_orders, assignedorders)
	for m in additionalorders
		sp_itemson[m] = copy(itemson[m])
	end

	sp_items = []
	for m in sp_orders
		sp_items = union(sp_items, sp_itemson[m])
	end

	#Update tabu windows
	if (length(tabulist) >= maxtabu) & (tabulist != [])
		popfirst!(tabulist)
	end
	push!(tabulist, sp_winid)

	return sp_winid, sp_orders, sp_window, sp_pods, sp_itemson, sp_items, tabulist
	
end

#-----------------------------------------------------------------------------------------------------#

function selectlearnthenoptimizesubproblem(currpartition, currsol, windows, windowscontaining, windowidlookup, windowsduring, windowsynergy, targetnumorders, targetnumpods, tabulist, lastoptimizeddifference, sp_iter)

	assignedorders_all, unassignedorders = findallunassignedorders(currpartition, currsol)
	objective, x, alpha, zeta = constructbestsubproblem(currpartition, currsol, features, featurenums, featureinfo, windows, windowscontaining, windowidlookup, windowsduring, windowsynergy, unassignedorders, tabulist, lastoptimizeddifference, sp_iter)
	sp_winid, sp_orders, sp_window, sp_pods, sp_itemson, sp_items, tabulist = parseoptimizationoutput(windows, currpartition, currsol, x, alpha, zeta, tabulist)

	return sp_winid, sp_orders, sp_window, sp_pods, sp_itemson, sp_items, tabulist, objective

end