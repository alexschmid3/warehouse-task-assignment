
function writepickdistrib(filename, globalsolution)

    df = DataFrame(pod = [], picktype = [], picks = [])
   
    for w in workstations, t in times
        uniquepods = unique([p for (m,i,p) in globalsolution.itempodpicklist[w,t]])
        itempicks, orderpicks = Dict(), Dict()
        for p in uniquepods
            itempicks[p] = 0
            orderpicks[p] = []
        end
        for (m,i,p) in globalsolution.itempodpicklist[w,t]
            itempicks[p] += 1
            push!(orderpicks[p], m)
        end
        for p in uniquepods
            push!(df, [p "Items picked" itempicks[p] ])
            push!(df, [p "Unique orders picked" length(unique(orderpicks[p]))])
        end
    end 

    CSV.write(filename, df)

end

function writedistancedistrib(filename, globalsolution)

    df = DataFrame(pod = [], itempicks = [], disttraveled = [])
   
    for w in workstations, t in times
        uniquepods = unique([p for (m,i,p) in globalsolution.itempodpicklist[w,t]])
        itempicks, orderpicks = Dict(), Dict()
        for p in uniquepods
            itempicks[p] = 0
        end
        for (m,i,p) in globalsolution.itempodpicklist[w,t]
            itempicks[p] += 1
        end
        for p in uniquepods
            disttraveled = manhattandist(podstorageloc[p],w)
            push!(df, [p itempicks[p] disttraveled])
        end
    end 

    CSV.write(filename, df)

end
