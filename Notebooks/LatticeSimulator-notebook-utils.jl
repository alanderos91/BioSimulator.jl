using Plots
gr(fmt = :png, size = (800, 600), dpi = 128)

##### constants #####

# for reproducibility
const CONFIGURATION_SEED = 5357
const SIMULATION_SEED = 1919

# feel free to modify these values for the notebook
MAX_TIME    = 10.0
MAX_POPSIZE = 10_000
MAX_STEPS   = 100

##### function definitions #####

function generate_random_2Dpoints(xlim, ylim, saturation)
    area = (xlim[2] - xlim[1] + 1) * (ylim[2] - ylim[1] + 1)
    N = Int(ceil(area * saturation))

    list = Tuple{Int,Int}[]

    while length(list) < N
        point = (rand(xlim[1]:xlim[2]), rand(ylim[1]:ylim[2]))
        if point ∉ list
            push!(list, point)
        end
    end

    points = zeros(Int64, 2, length(list))
    for i in eachindex(list)
        points[1, i] = list[i][1]
        points[2, i] = list[i][2]
    end

    return points
end

function plot_config(lattice; show_open = false, title = "", xlab = "", ylab = "", labels = nothing)
    # get number of particle types
    L = lattice.L
    
    # set labels if they were omitted
    if labels == nothing
        labels = ["type $(i)" for i in 1:L]
    end
    
    # retrieve the coordinates of every site tracked by the lattice
    sites = collect(values(lattice.coord2site))

    # retrieve the types of every site
    point_types = Vector{Vector{Tuple{Int,Int}}}(L + 1)

    # extract open sites
    point_types[1] = [tuple(site...) for site in filter(y -> get_ptype(y) == 1, sites)]
    
    # extract remaining sites
    for i in 2:L+1
        point_types[i] = [tuple(site...) for site in filter(y -> get_ptype(y) == i, sites)]
    end
    
    # initialize the plot
    p = plot(xlab = xlab, ylab = ylab, title = title)

    # should we display the open sites?
    if show_open
        scatter!(point_types[1], label = "open", color = :white)
    end
    
    # add the normal types
    for i in 2:L+1
        scatter!(point_types[i], label = labels[i - 1])
    end

    return p 
end



function plot_hex_config(lattice; show_open = false, title = "", xlab = "", ylab = "", labels = nothing)
    # get number of particle types
    L = lattice.L
    
    # set labels if they were omitted
    if labels == nothing
        labels = ["type $(i)" for i in 1:L]
    end
    
    # retrieve the coordinates of every site tracked by the lattice
    sites = collect(values(lattice.coord2site))
    
    new_coords = zeros(Float64, 2, length(sites))

    thingx = sqrt(1 / 4 - 1 / 16)
    thingy = 3 / 4

    for i in 1:length(sites)
        # janky version, modify the sites by shifting the (x,y) coordinates over by the appropriate amount.  Depends on trig.   Done in axial coordinates
        x = sites[i].coord[1]
        y = sites[i].coord[2]
        #x = coordinates(sites[i])[1]
        #y = coordinates(sites[i])[2]

        new_coords[1, i] = x + thingx * y
        new_coords[2, i] = thingy * y

    end



    # retrieve the types of every site
    #point_types = Vector{Vector{Tuple{Float64,Float64}}}(L + 1)


    # extract open sites
    #point_types[1] = [(site.coord[1] + thingx * site.coord[2], thingy * site.coord[2]) for site in filter(y -> get_ptype(y) == 1, sites)]
    
    # extract remaining sites
   # for i in 2:L+1
   #     point_types[i] = [(coordinates(site)[1] + thingx * coordinates(site)[2], thingy * coordinates(site)[2]) for site in filter(y -> get_ptype(y) == i, sites)]
   # end
    
    # initialize the plot
    p = plot(xlab = xlab, ylab = ylab, title = title)

    # should we display the open sites?
   # if show_open
   #     scatter!(point_types[1], label = "open", color = :white)
   # end
    
    # add the normal types
    #for i in 2:L+1
    for i in 2:L+1
            scatter!( [(site.coord[1] + thingx * site.coord[2], thingy * site.coord[2]) for site in filter(y -> y.state.ptype == i, sites)], label = labels[i - 1], 
            markershape = :hexagon, markersize = 2.5, markerstrokewidth = 0.3)
    end

    return p 
end




function generate_random_2Dhex(xlim, ylim, saturation)
    area = (xlim[2] - xlim[1] + 1) * (ylim[2] - ylim[1] + 1)
    N = Int(ceil(area * saturation))

    list = Tuple{Int,Int}[]

    while length(list) < N
        point = (rand(xlim[1]:xlim[2]), rand(ylim[1]:ylim[2]))
        if point ∉ list
            push!(list, point)
        end
    end

    points = zeros(Int64, 2, length(list))
    for i in eachindex(list)
        points[1, i] = list[i][1]
        points[2, i] = list[i][2]
    end

    return points
end

function generate_bordered_2D_hexgrid(l, saturation)
    # l is the radius of the hexagon, l >= 2
    
    # number of tiles in a hexagon of length l
    area = 3 * l * (l + 1) + 1

    N = Int(ceil(area * saturation)) - (6 * l) # take off the border

    list = Tuple{Int,Int}[]

    # add in the border cells first; there are 6*l of these

    for x in 0:(l - 1)
        point = (- l , x)
        push!(list, point)

        point = (- l + x , l)
        push!(list, point)

        point = (x , l - x) 
        push!(list, point)

        point = (l , - x)
        push!(list, point)

        point = (l - x , - l)
        push!(list, point)

        point = (- x , - l + x)
        push!(list, point)

    end

    # while length(list) < (N + 6 * l)
    #     point = (rand(xlim[1]:xlim[2]), rand(ylim[1]:ylim[2]))
    #     if point ∉ list
    #         push!(list, point)
    #     end
    # end

    points = zeros(Int64, 2, length(list))
    for i in eachindex(list)
        points[1, i] = list[i][1]
        points[2, i] = list[i][2]
    end

    return points
end


