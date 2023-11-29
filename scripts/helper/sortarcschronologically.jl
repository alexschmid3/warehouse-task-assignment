function sortarcschronologically(arclist)
    return sort(arclist, by=x->nodelookup[arclookup[x][1]][2])
end