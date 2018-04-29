struct PetriNet
  g :: DiGraph

  species_nodes  :: Vector{Int}
  reaction_nodes :: Vector{Int}

  species_styles  :: Dict{Int,String}
  reaction_styles :: Dict{Int,String}

  node_labels :: Vector{String}
  edge_labels :: Dict{Tuple{Int,Int},String}
  edge_styles :: Dict{Tuple{Int,Int},String}
end

function petri_net(model :: Network)
  s = n_species(model)
  r = n_reactions(model)
  g = DiGraph(s + r)

  species   = species_list(model)
  reactions = reaction_list(model)
  id2ind    = Dict( id => i for (i, id) in enumerate(keys(species)) )
  edge_set = Tuple{Int,Int}[]
  stoc_set = Dict{Tuple{Int,Int},Int}()

  edge_styles = Dict{Tuple{Int,Int},String}()

  # construct edges and extract stoichiometries
  # species numbered 1 thru s
  # reactions numbered s+1 thru s+r
  for (k, reaction) in enumerate(values(reactions))
    j = k + s
    reactants = reaction.reactants
    products  = reaction.products

    for (reactant, v) in reactants
      i = id2ind[reactant]
      e = (i, j)
      push!(edge_set, e)
      edge_styles[e] = "-stealth, draw, rounded corners=5pt, solid, xshift=-2pt"
      v > 1 && (stoc_set[e] = v)
    end

    for (product, v) in products
      i = id2ind[product]
      e = (j, i)
      push!(edge_set, e)
      edge_styles[e] = "-stealth, draw, rounded corners=5pt, dashed, xshift=2pt, red, swap"
      v > 1 && (stoc_set[e] = v)
    end
  end

  for edge in edge_set
    add_edge!(g, edge[1], edge[2])
  end

  species_nodes  = collect(1:s)
  species_labels = map(string, keys(species))
  species_styles = Dict( i => "draw, rounded corners, fill=blue!10" for i in species_nodes )

  reaction_nodes = collect(s+1:s+r)
  reaction_labels = map(string, keys(reactions))
  reaction_styles = Dict( i => "draw, rounded corners, thick, fill=red!10" for i in reaction_nodes)

  node_labels = [ species_labels; reaction_labels ]
  edge_labels = Dict( e => string(v) for (e, v) in stoc_set )

  return PetriNet(
    g,
    species_nodes,
    reaction_nodes,
    species_styles,
    reaction_styles,
    node_labels,
    edge_labels,
    edge_styles
  )
end

function draw(x :: PetriNet, options)
  graph           = x.g
  labels          = x.node_labels
  species_styles  = x.species_styles
  reaction_styles = x.reaction_styles
  edge_labels     = x.edge_labels
  edge_styles     = x.edge_styles

  labels = map(x -> replace(x, "_", "\$\\cdot\$"), labels)

  TikzGraphs.plot(graph, TikzGraphs.Layouts.Layered(),
    labels,
    node_styles = merge(species_styles, reaction_styles),
    edge_labels = edge_labels,
    edge_styles = edge_styles,
    #edge_style = "out=20, relative=true",
    options = options
  )
end

visualize(x :: Network, options = "") = draw(petri_net(x), options)
