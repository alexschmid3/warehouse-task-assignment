

using CSV, DataFrames, Graphs, GraphPlot, Compose

filename = "outputs/ordergraph/run207_2024-07-20/ordergraph.csv"

#Read data
data = CSV.read(filename, DataFrame)
workstations = unique(data[:,1])
times = unique(data[:,2])
orders = unique(data[:,3])
pods = unique(data[:,4])
items = unique(data[:,5])

#Map orders to indices
mapordertoindex, mapindextoorder = Dict(), Dict()
orderindex = 1
for m in orders
    mapordertoindex[m] = orderindex
    mapindextoorder[orderindex] = m
    orderindex += 1
end
numorders = length(orders)

#Initialize graph
g=SimpleGraph(numorders,0)

#Add edges between orders that share a pod
podpicks, orderspulledfrom = Dict(), Dict()
for row in 1:size(data)[1]
    w,t,m,p,i = data[row,1], data[row,2], data[row,3], data[row,4], data[row,5]
    orderspulledfrom[p,w,t] = []
    push!(podpicks, (p,w,t))
end
for row in 1:size(data)[1]
    w,t,m,p,i = data[row,1], data[row,2], data[row,3], data[row,4], data[row,5]
    push!(orderspulledfrom[p,w,t], m)
end
for (p,w,t) in podpicks
    for m1 in orderspulledfrom[p,w,t], m2 in setdiff(orderspulledfrom[p,w,t], m1)
        add_edge!(g, mapordertoindex[m1], mapordertoindex[m2])
    end
end

add_edge!(g, 1, 6);


draw(PNG("tempgraph.png", 16cm, 16cm), gplot(g))
