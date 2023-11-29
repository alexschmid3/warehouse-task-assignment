
function updatesolvemetrics(s, solvetime, solvetime_sp, spselectiontime)

	solvemetrics.solve_time[s] += solvetime
	solvemetrics.solvetime_sp[s] += solvetime_sp
	solvemetrics.lsnsiterations[s] += 1
	solvemetrics.solvetime_spsel[s] += spselectiontime

end