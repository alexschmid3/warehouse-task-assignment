function getcurrentobjective(currpartition, currsol)
    obj = sum(sum(length(currsol.itempodpicklist[w,t]) for t in times) for w in currpartition.workstations)
    println("Throughput = ", obj)
end