function parsetabucodes(tabutype)

	if tabutype == "orderedsp"
		windowforcingflag = 1
		maxtabu = 0
		lastoptpenaltyflag = 0
	elseif tabutype == "tabu8"
		windowforcingflag = 0
		maxtabu = 8
		lastoptpenaltyflag = 0
	elseif tabutype == "tabu15"
		windowforcingflag = 0
		maxtabu = 15
		lastoptpenaltyflag = 0
	elseif tabutype == "lastopt10"
		windowforcingflag = 0
		maxtabu = 10
		lastoptpenaltyflag = 1
	elseif tabutype == "lastoptall"
		windowforcingflag = 0
		maxtabu = 100
		lastoptpenaltyflag = 1
	end

	return windowforcingflag, maxtabu, lastoptpenaltyflag

end

