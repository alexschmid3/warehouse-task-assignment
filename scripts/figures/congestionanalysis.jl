
function congestionanalysis(nocongestion_flag)

    l = maps.mapintersectiontorow[65]
    totalcong = sum(currcong[p] for p in  partitioninfo[1].pods)

    allcongestionutils = []
    intersectionviolations, maxtraffic = Dict(), Dict()
    for i in globalpartition.intersections[1]
        intersectionviolations[i] = 0
        maxtraffic[i] = 0
    end

    for s in 1:numpartitions
        currpartition = partitioninfo[s]
        for i in currpartition.intersections, t in 0:congestiontstep:horizon
            totalutil = sum(currcong[p][maps.mapintersectiontorow[i],maps.maptimetocolumn[t]] for p in currpartition.pods) / intersectionmaxpods[i]
            push!(allcongestionutils, totalutil)
            maxtraffic[i] = max(maxtraffic[i], totalutil)
            if totalutil > 1 + 1e-4
                intersectionviolations[i] += 1
            end
        end
    end
    
    df = DataFrame(congutil = allcongestionutils)

    if nocongestion_flag == 0
        CSV.write("figures/congestionanalysis/congestionsolution2.csv", df)
    else nocongestion_flag == 1
        CSV.write("figures/congestionanalysis/nocongestionsolution2.csv", df)
    end

    include("scripts/visualizations/warehouseviz.jl")
    warehouseviz("figures/congestionanalysis/warehousecongestion_nocong.png", 4000)

end
