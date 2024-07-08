function formatoutcomes(oldobj_sp, newobj_sp, solvetime_sp)

	improvement = newobj_sp - oldobj_sp
	if (oldobj_sp == 0) & (newobj_sp == 0)
		improvementpct = 0
	elseif (oldobj_sp == 0)
		improvementpct = 1
	else
		improvementpct = (newobj_sp - oldobj_sp)/oldobj_sp
	end
	
	outcomes = [oldobj_sp newobj_sp improvement improvementpct solvetime_sp]

	return outcomes

end