
function parsemethodname(methodname)
	
	if occursin("_",methodname)
		shortmethodname = methodname[1:findfirst("_",methodname)[1]-1]
	else
		shortmethodname = methodname
	end
	if shortmethodname == "learnbench"
		subproblembudget = parse(Int,methodname[findfirst("_",methodname)[1]+1:length(methodname)])
	else
		subproblembudget = 0
	end

	return shortmethodname, subproblembudget

end
