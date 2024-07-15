
function congestionanalysis(nocongestion_flag)

    allcongestionutils = []

    for s in 1:numpartitions
        currpartition = partitioninfo[s]
        for i in currpartition.intersections, t in 0:congestiontstep:horizon
            push!(allcongestionutils, sum(currcong[p][maps.mapintersectiontorow[i],maps.maptimetocolumn[t]] for p in currpartition.pods) / intersectionmaxpods[i])
        end
    end
    
    df = DataFrame(congutil = allcongestionutils)
    
    if nocongestion_flag == 0
        CSV.write("figures/congestionanalysis/congestionsolution.csv", df)
    else nocongestion_flag == 1
        CSV.write("figures/congestionanalysis/nocongestionsolution.csv", df)
    end

end