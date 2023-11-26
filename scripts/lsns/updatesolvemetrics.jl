
function updatesolvemetrics(s, solvetime, solvetime_sp)

	solvemetrics.solve_time[s] += solvetime
	solvemetrics.solvetime_sp[s] += solvetime_sp
	solvemetrics.lsnsiterations[s] += 1

end