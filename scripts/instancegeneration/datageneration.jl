using JuMP
using Gurobi
using CSV
using DataFramesMeta
using Statistics
using StatsPlots
using MLBase
using GLM
using Distances
using Plots.PlotMeasures
using PlotlyBase
using ParallelKMeans
using Random
using Distributions
using Distances
using DataFrames;

##Input:
#xm = warehouse length in x dimension (meters)
#ym = warehouse length in y dimension (meters)
#Workstations = # of workstations

#Px = pods per block in x direction
#Py = pods per block in y direction
#Street_Width = width of streets

function Warehouse(xm,ym,Workstations,Seed;Px=6,Py=4,Street_Width=2)
    
    #xm,ym,Workstations,Seed =  warehouse_x_length_meters,  warehouse_y_length_meters, num_workstations, random_seed

    Group_width_x = Px
    Group_width_y = Py

    #x and y length of a group/street pair:
    Street_Group_x = Group_width_x + Street_Width
    Street_Group_y = Group_width_y + Street_Width

    #Rows in x and y direction
    x_rows = convert(Int64, floor(xm/Street_Group_x))
    y_rows = convert(Int64, floor(ym/Street_Group_y))
    if ym > 45
        storageloc_rows = 2:y_rows-1
    else
        storageloc_rows = 1:y_rows-1
    end

    #For loops to define streets:
    Vert_Streets = DataFrame(Array{Float64}(undef, x_rows, 2),[:left,:right])
    Hor_Streets = DataFrame(Array{Float64}(undef, y_rows, 2),[:top,:bottom])

    for i = 1:x_rows
        Vert_Streets[i,1] = (i-1) * Street_Group_x
        Vert_Streets[i,2] = (i-1) * Street_Group_x + Street_Width
    end
    for j = 1:y_rows
        Hor_Streets[j,1] = (j-1) * Street_Group_y
        Hor_Streets[j,2] = (j-1) * Street_Group_y + Street_Width
    end

    #Pod "Group" Locations (Bottom Right Corner):
    PodLocations = DataFrame(Array{Float64}(undef, x_rows * (y_rows - 1), 2),[:x,:y])
    count=1
    for i = 1:x_rows
        for j in storageloc_rows
            PodLocations[count,1] = i*Street_Width + (i-1)*Px
            PodLocations[count,2] = j*Street_Width + (j-1)*Py
            count=count+1
        end
    end

    #Workstation Locations 
    Work_Loc = DataFrame(Array{Float64}(undef, Workstations, 3),[:Station,:x,:y])
    Work_Loc.Station = collect(size(PodLocations)[1]+1:size(PodLocations)[1]+Workstations)

    #Big warehouse (bottom row and top row, equally spread)
    if ym > 45
        topStations, bottomStations = 1:convert(Int,ceil(Workstations/2)), convert(Int, ceil(Workstations/2)+1):Workstations
        stationsperpodcolumn = floor(x_rows / ceil(Workstations/2))
        for i in topStations
            Work_Loc[i,2] = Street_Width + (i-1)*(Street_Width + Group_width_x)*stationsperpodcolumn
            Work_Loc[i,3] = 1*Street_Width + (1-1)*Py
        end
        for i in bottomStations
            Work_Loc[i,2] = Street_Width + (i-length(topStations)-1)*(Street_Width + Group_width_x)*stationsperpodcolumn
            Work_Loc[i,3] = y_rows*Street_Width + (y_rows-1)*Py
        end
    #Small warehouse (bottom row, equally spread)
    else
        bottomStations = 1:Workstations
        stationsperpodcolumn = floor(x_rows / Workstations)
        for i in bottomStations
            Work_Loc[i,2] = Street_Width + (i-1)*(Street_Width + Group_width_x)*stationsperpodcolumn
            Work_Loc[i,3] = y_rows*Street_Width + (y_rows-1)*Py
        end
    end

    #Number of Pods We can fit:
    Num_Pods = (x_rows * length(storageloc_rows)) * Px * Py

    #Pod location dataframe
    Pod_Start = DataFrame(Array{Float64}(undef, Num_Pods, 3),[:Pod,:start_x,:start_y])
    Pod_Start.Pod = collect(1:Num_Pods)
    for i = 1:(x_rows * length(storageloc_rows))
        Pod_Start[((i-1)*(Px * Py)+1):(i*(Px * Py)),2] = repeat([PodLocations[i,1]],Px * Py)
        Pod_Start[((i-1)*(Px * Py)+1):(i*(Px * Py)),3] = repeat([PodLocations[i,2]],Px * Py)
    end

    return(Vert_Streets, Hor_Streets, Work_Loc, Pod_Start)

end

#N = # of unique items in the warehouse (Going to set this at 1/16 of capacity)
#n = # of items held by each pod (We hold this constant for now, guess: 200 (~33 shelves per pod))
#Pods = # of pods in our problem

function Inventory(N, n, Pods, Seed)
    
    #Normal distribution centered at N/2 where 0 is 3 sd out
    Random.seed!(Seed+N)
    x = rand(Normal(N/2, (N - N/2)/3), n*Pods-N)
    x = convert.(Int64, round.(x))
    
    #Items = DataFrame(Array{Float64}(undef, N, 2),[:Item, :Inventory])
    #Items.Item = collect(1:N)
    Items = DataFrame(Array{Float64}(undef, n*Pods, 2),[:Item, :Inventory])
    Items.Inventory = repeat([0.0],n*Pods)
    
    #Assume there is inventory of one of at least everything
    for i = 1:N
            Items.Item[i] = i
            Items.Inventory[i] = 1
    end
    
    #The rest is based on the normal distribution draws.
    for i = 1: size(x)[1]
        if 0 < x[i] && x[i] <= N
            Items.Item[i+N] = x[i]
            Items.Inventory[i+N] = 1
        end
    end
    
    Items = Items[(Items[!, :Inventory].!= 0), :]
    
    #Now our inventory dataframe    
    Inventory = DataFrame(Array{Float64}(undef, 0, 3),[:Pod, :Item, :Inventory])
    
    Random.seed!(Seed+N)
    order = shuffle(collect(1:size(Items)[1]))
    
    #Best
    for i = 1:Pods
        for j = 1:n
            if ((i-1)*n + j) <= size(Items)[1]
                Item = Items[order[(i-1)*n + j],1]
                push!(Inventory,(i,Item,1))
            end
        end
    end  
    Inventory = DataFrames.groupby(Inventory, [:Pod, :Item])
    Inventory = DataFrames.combine(Inventory, :Inventory => sum)

    return(Inventory)
    
end

#per_hr = Orders per hour (at full scale, this is over 1,000)
#p = probability of one item orders
#q = Geometric distribution parameter for multi-item orders
#N = # of unique items in the warehouse (Going to set this at 1/16 of capacity)

#T = full time window (in seconds)
#comp_time = Seconds per item to complete the order

function Orders(per_hr, p, q, N, Seed; comp_time = 1800, T=43200)
    
    #The random seed gets set a bunch of times in this script because for some reason, this is the only way to replicate the instance correctly on the cluster
    #Otherwise, you end up with different orders when you run locally vs. remotely
    #Not sure why, but thanks for figuring that out, Riley!

    Random.seed!(2*Seed)

    num_orders = (per_hr/3600) * T
    num_orders = convert(Int64, floor(num_orders))
    
    Order_Details = DataFrame(Array{Float64}(undef, num_orders, 4),[:Order,:num_items,:arrival,:due])
    Order_Details.Order = collect(1:num_orders)
    
    #Finding # of items per order:
    x = repeat([0], num_orders)
    
    Random.seed!(Seed + num_orders)
    y = rand(Uniform(0,1), num_orders)

    counter1, counter2 = 0, 0 
    for i = 1:num_orders
        if y[i] <= p
            x[i] = 1
            counter1 += 1
        else
            Random.seed!(Seed + i)
            x[i] = rand(Geometric(q), 1)[1] .+ 2
        end
        counter2 += 1
    end
    println("Percent = ", counter1/ counter2)

    Order_Details.num_items = x
    
    #Arrival and Due times
    Random.seed!(3*Seed)
    Order_Details.arrival = sample(1:T,num_orders)
    Order_Details.due = Order_Details.arrival + comp_time .* Order_Details.num_items
    
    #Now starting on itemized order list
    Itemized = DataFrame(Array{Float64}(undef, convert(Int64,sum(Order_Details.num_items)), 2),[:Order,:Item])
    Random.seed!(Seed+N)
    y = rand(Normal(N/2, (N - N/2)/3), convert(Int64, round(1.3 * sum(Order_Details.num_items))))
    y = convert.(Int64, round.(y))
    y = y[0 .< y]
    y = y[y .<= N]
    y = y[1:convert(Int64,sum(Order_Details.num_items)),:]
    
    #Formatting the data frame
    count = 0
    for i = 1:num_orders
        Itemized.Order[(count+1):(count+1convert(Int64,Order_Details.num_items[i]))] = repeat([Order_Details.Order[i]],convert(Int64,Order_Details.num_items[i]))
        count += convert(Int64,Order_Details.num_items[i])
    end
    for i = 1:size(Itemized,1)
        Itemized.Item[i] = y[i]
    end
    
    return (Order_Details,Itemized)
    
end

##Input:
#xm = warehouse length in x dimension (meters)
#ym = warehouse length in y dimension (meters)
#Workstations = # of workstations
#N = # of unique items in the warehouse (Going to set this at 1/16 of capacity)
#n = # of items held by each pod (We hold this constant for now, guess: 200 (~33 shelves per pod))
#Pods = # of pods in our problem
#per_hr = Orders per hour (at full scale, this is over 1,000)
#p = probability of one item orders
#q = Geometric distribution parameter for multi-item orders

#Px = pods per block in x direction
#Py = pods per block in y direction
#Street_Width = width of streets
#T = full time window (in seconds)
#comp_time = Seconds allowed per item to complete the order

#xm, ym, Workstations, N, n, per_hr, p, q, Seed = warehouse_x_length_meters, warehouse_y_length_meters, num_workstations, num_unique_items, num_items_per_pod, num_orders_per_hour, prob_one_item_orders, geo_dist_param_orders, random_seed
#Px, Py,Street_Width,comp_time, T=6, 4, 2, 1800, 43200

function Data(xm, ym, Workstations, N, n, per_hr, p, q, Seed;Px=6,Py=4,Street_Width=2,comp_time = 1800, T=43200)
    
    Random.seed!(Seed)

    Ware = Warehouse(xm,ym,Workstations,Seed;Px,Py,Street_Width)
    Inv = Inventory(N,n,size(Ware[4],1),Seed)
    Ord = Orders(per_hr,p,q,N,Seed;comp_time,T)

    #Re-number orders
    
    
    #Output Order:
    #1) Vertical Streets
    #2) Horizontal Streets
    #3) Workstation Locations
    #4) Pod "Homes"
    #5) Inventory
    #7) Order Details
    #8) Itemized Order List
    
    return(Ware,Inv,Ord)
    
end

#xm, ym, Workstations, N, n, per_hr, p, q, Seed = warehouse_x_length_meters, warehouse_y_length_meters, num_workstations, num_unique_items, num_items_per_pod, num_orders_per_hour, prob_one_item_orders, geo_dist_param_orders, random_seed