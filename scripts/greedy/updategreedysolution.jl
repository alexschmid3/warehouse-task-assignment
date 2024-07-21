
#-----------------------------------------------------------------------------#

function update_hv_greedy(currsol, currsol_greedysets, gp, h)

	currsol.stationassign[gp.m] = gp.w
    
    earliestt, latestt, orderclosed = min(horizon, gp.orderendtime), max(0, gp.orderstarttime), 1
    for i in itemson[gp.m]
        mysum = 0
        for p in intersect(gp.relevantpods, currsol_greedysets.podswith_greedy[i]), t in max(0,gp.podstarttime[p]):tstep:min(horizon,gp.podendtime[p])
            if value(h[i,p,t]) > 1e-4
                currsol_greedysets.remaininginventory[i, p] -= 1
                earliestt = min(earliestt, t)
                latestt = max(latestt, t) 
                currsol.podsworkedat[gp.w,t] = union(currsol.podsworkedat[gp.w,t], p)
                remove!(currsol_greedysets.podswith_greedy[i], p)
                push!(currsol.itempodpicklist[gp.w, t], (gp.m, i, p))
                for t_prime in t:tstep:horizon
                    currsol.h[gp.m, i, p, gp.w, t_prime] = 1
                end
                mysum += value(h[i,p,t])
                break
            end
        end
        #if sum(sum(value(h[i,p,t]) for t in gp.podstarttime[p]:tstep:gp.podendtime[p]; init=0) for p in intersect(gp.relevantpods, currsol_greedysets.podswith_greedy[i]); init=0) < 1e-4  
        #for p in intersect(gp.relevantpods, currsol_greedysets.podswith_greedy[i]), t in gp.podstarttime[p]:tstep:gp.podendtime[p]
        #    mysum += value(h[i,p,t])
        #end
        if mysum < 1e-4
            orderclosed = 0
            latestt = horizon
        end
    end

    #Update the open orders with the assignment of m
    #Remember, we don't consider an order open in the period it is closed
    for t in earliestt:tstep:latestt-tstep*orderclosed
        push!(currsol.ordersopen[gp.w, t], gp.m)
        currsol.v[gp.m, gp.w, t] = 1
    end

	return currsol, currsol_greedysets

end

#-----------------------------------------------------------------------------#

function update_y_greedy(currsol, gp, y, z)

    for p in gp.relevantpods
        podtostationtime = arclength[podstorageloc[p], gp.w]
        awaystart, awayend, podarrives = 0, 0, 0
        for t in gp.podstarttime[p]:tstep:gp.podendtime[p]
            if value(y[p,t]) > 0.01
                awaystart += max(-1*tstep, t - podtostationtime)
                podarrives += t
            end
            if value(z[p,t]) > 0.01
                awayend += min(horizon, t + podtostationtime)
                for t_p in awaystart:tstep:awayend
                    push!(currsol.podbusy[p], (gp.w, t_p))
                end

                #Pod trip to the station
                if awaystart < 0
                    a = extendedarcs[extendednodes[podstorageloc[p], dummystarttime], extendednodes[gp.w, podarrives]]
                    currsol.y[p, a] = 1
                    push!(currsol.ypath[p], a)
                else
                    a = arcs[nodes[podstorageloc[p], awaystart], nodes[gp.w, awaystart + podtostationtime]]
                    currsol.y[p, a] = 1
                    push!(currsol.ypath[p], a)
                end		
                                    
                #Pod trip back to storage region
                if t + podtostationtime > horizon
                    a = extendedarcs[extendednodes[gp.w, t], extendednodes[podstorageloc[p], dummyendtime]]
                    currsol.y[p, a] = 1	
                    push!(currsol.ypath[p], a)
                else 
                    a = arcs[nodes[gp.w, awayend - podtostationtime], nodes[podstorageloc[p], awayend]]
                    currsol.y[p, a] = 1	
                    push!(currsol.ypath[p], a)
                end	

                #Pod stay in the workstation queue
                for t_p in awaystart+podtostationtime:tstep:min(horizon-tstep, t-tstep)
                    a = arcs[nodes[gp.w, t_p], nodes[gp.w, t_p + tstep]]
                    currsol.y[p,a] = 1
                    push!(currsol.ypath[p], a)
                end

                #Pod no longer stationed at home
                for t_p in awaystart:tstep:min(dummyendtime-tstep, awayend - tstep)
                    a = extendedarcs[extendednodes[podstorageloc[p], t_p], extendednodes[podstorageloc[p], t_p + tstep]]
                    currsol.y[p, a] = 0
                    remove!(currsol.ypath[p], a)
                end

                currsol.ypath[p] = sortarcschronologically(currsol.ypath[p])
                break
            end
        end
    end

	return currsol

end

#-----------------------------------------------------------------------------#

function update_picklists_greedy(currsol, currsol_greedysets, gp, h)

	#Update item-pod pick list and podsworked based on the new picks
    for i in itemson[gp.m]
        for p in intersect(gp.relevantpods, currsol_greedysets.podswith_greedy[i]), t in gp.podstarttime[p]:tstep:gp.podendtime[p]
            if value(h[i,p,t]) > 0.01
                push!(currsol.itempodpicklist[gp.w,t], (gp.m,i,p))
                currsol.podsworkedat[gp.w,t] = union(currsol.podsworkedat[gp.w,t], p)
                break
            end
        end
    end

	return currsol

end

#-----------------------------------------------------------------------------#

function updateorders_greedy(currsol, gp)

	remove!(currsol.unassignedorders, gp.m)

	return currsol

end

#-----------------------------------------------------------------------------#

function updatecongestion_greedy(currsol, gp)

    for p in gp.relevantpods
        currcong[p] = sum(congestionsignature[a] for a in [a2 for a2 in currsol.ypath[p] if a2 <= numarcs])
    end

end

#-----------------------------------------------------------------------------#

#currsol, currsol_greedysets, gp, currpartition, h, y, z = currsol, currsol_greedysets, gp, currpartition, h_assign, y_assign, z_assign 

function updategreedysolution(currsol, currsol_greedysets, gp, currpartition, h, y, z)

    currsol, currsol_greedysets = update_hv_greedy(currsol, currsol_greedysets, gp, h)
    currsol = update_y_greedy(currsol, gp, y, z)
    currsol = update_picklists_greedy(currsol, currsol_greedysets, gp, h)
    currsol = updateorders_greedy(currsol, gp)
    updatecongestion_greedy(currsol, gp)

	return currsol, currsol_greedysets

end

