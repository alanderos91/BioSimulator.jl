struct NeighborhoodClass{S}
  members::Vector{S}
  class::Int
end

function NeighborhoodClass{S}(class) where S
  members = S[]; sizehint!(members, 100)

  return NeighborhoodClass{S}(members, class)
end

# iteration #
iterate(nc::NeighborhoodClass) = iterate(nc.members)
iterate(nc::NeighborhoodClass, state) = iterate(nc.members, state)

# indexing #

# start(nc::NeighborhoodClass) = start(nc.members)
# next(nc::NeighborhoodClass, state) = next(nc.members, state)
# done(nc::NeighborhoodClass, state) = done(nc.members, state)

getindex(nc::NeighborhoodClass, i) = getindex(nc.members, i)
setindex!(nc::NeighborhoodClass, val, i) = setindex!(nc.members, val, i)
eachindex(nc::NeighborhoodClass) = eachindex(nc.members)

# sampling #

function sample_center(nc::NeighborhoodClass, l, nl)
  x = nc.members
  c = nl * rand()
  j = 1
  s = get_ptype(x[j]) == l ? 1 : 0
  
  while s < c && j < length(x)
    isequal(get_ptype(x[j += 1]), l) && (s += 1)
  end

  return x[j]
end

function add_member!(nc::NeighborhoodClass, x::Site)
  if nc.class == get_neighbor_class(x)
    # i = searchsortedfirst(nc.members, x, by = label)
    i = searchsortedfirst(nc.members, x) # for v0.6
    insert!(nc.members, i, x)
  else
    error("class of $(x) != $(nc.class)")
  end

  return nothing
end

function rmv_member!(nc::NeighborhoodClass, x::Site)
  # i = searchsortedfirst(nc.members, x, by = label)
  i = searchsortedfirst(nc.members, x) # for v0.6
  deleteat!(nc.members, i)

  return nothing
end

function Base.show(io::IO, nc::NeighborhoodClass)
  print(io, "----- Class ", nc.class, "-----\n")
  print(io, nc.members)
  return nothing
end
