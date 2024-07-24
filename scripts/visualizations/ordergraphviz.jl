
using CSV, DataFrames, Graphs, Plots, GraphRecipes, Cairo, Fontconfig

#Data
#filename = "outputs/ordergraph/run203_2024-07-20/ordergraph.csv" #LTO
filename = "outputs/ordergraph/run207_2024-07-20/ordergraph.csv" #Synergy - Alone
#filename = "outputs/ordergraph/run208_2024-07-21/ordergraph.csv" #Synergy - Greedy
#filename = "outputs/ordergraph/run209_2024-07-21/ordergraph.csv" #Random
#filename = "outputs/ordergraph/run210_2024-07-21/ordergraph.csv" #Greedy

#Graph design
colorlist = [colorant"#648FFF",colorant"#FE6100",colorant"#DC267F"]

#Read data
data = CSV.read(filename, DataFrame)
workstations = unique(data[:,1])
times = unique(data[:,2])
orders = unique(data[:,3])
pods = unique(data[:,4])
items = unique(data[:,5])

workstationcolor, workstationcolor_halo = Dict(), Dict()
for w in workstations
    workstationcolor[w] = colorlist[w - minimum(workstations) + 1]
    #workstationcolor_halo[w] = halocolorlist[w - minimum(workstations) + 1]
end

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
podpicks, orderspulledfrom = [], Dict() 
nodelabels = [" " for m in 1:numorders]
halos = zeros(numorders)
ordercolors = Array{Any}(undef, numorders, 1)
ordersize = zeros(numorders)
for row in 1:size(data)[1]
    w,t,m,p,i = data[row,1], data[row,2], data[row,3], data[row,4], data[row,5]
    orderspulledfrom[(p,w,t)] = []
    push!(podpicks, (p,w,t))
    ordercolors[mapordertoindex[m]] = workstationcolor[w]
    ordersize[mapordertoindex[m]] += 1
end
for row in 1:size(data)[1]
    w,t,m,p,i = data[row,1], data[row,2], data[row,3], data[row,4], data[row,5]
    push!(orderspulledfrom[p,w,t], m)
end
for (p,w,t) in unique(podpicks)
    for m1 in orderspulledfrom[p,w,t], m2 in setdiff(orderspulledfrom[p,w,t], m1)
        add_edge!(g, mapordertoindex[m1], mapordertoindex[m2])
    end
    for m1 in unique(orderspulledfrom[p,w,t])
        if length([m for m in orderspulledfrom[p,w,t] if m==m1]) > 1.001
            println("orderspulledfrom[$p,$w,$t] = ", orderspulledfrom[p,w,t])
            halos[mapordertoindex[m1]] += length([m for m in orderspulledfrom[p,w,t] if m==m1]) - 1
            nodelabels[mapordertoindex[m1]] = string(length([m for m in orderspulledfrom[p,w,t] if m==m1]) - 1) #"+"
            #add_edge!(g, mapordertoindex[m1], mapordertoindex[m1])
        end
    end
end

ordersize = ordersize./10

println("Clustering coefficient = ", global_clustering_coefficient(g))
println("Single-order consolidation = ", sum(halos))
println("Cross-order consolidation = ", ne(g))

#Create graph
graphplot(g,
          nodeshape=:circle, 
          nodesize=0.15,
          axis_buffer=0.1,
          curves=false,
          color=:black,
          nodecolor=ordercolors,
          names=nodelabels,
          fontsize=7,
          linewidth=1)