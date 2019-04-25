struct SampleClass{S}
  members::Vector{S}
  class::Int
end

function SampleClass{S}(class) where S
  members = S[]; sizehint!(members, 100)

  return SampleClass{S}(members, class)
end

# iteration #

iterate(nc::SampleClass) = iterate(nc.members)
iterate(nc::SampleClass, state) = iterate(nc.members, state)

# indexing #

getindex(nc::SampleClass, i) = getindex(nc.members, i)
setindex!(nc::SampleClass, val, i) = setindex!(nc.members, val, i)
eachindex(nc::SampleClass) = eachindex(nc.members)

# sampling #

# overload StatsBase.sample
sample(nc::SampleClass) = sample(nc.members)

function add_member!(nc::SampleClass, x::Site)
  # if nc.class == get_neighbor_class(x)
  #   push!(nc.members, x)
  # else
  #   error("class of $(x) != $(nc.class)")
  # end
  # i = searchsortedfirst(nc.members, x, by = label)
  i = searchsortedfirst(nc.members, x) # for v0.6
  insert!(nc.members, i, x)
  
  return nothing
end

function rmv_member!(nc::SampleClass, x::Site)
  # i = searchsortedfirst(nc.members, x, by = label)
  i = searchsortedfirst(nc.members, x) # for v0.6
  deleteat!(nc.members, i)

  return nothing
end

function Base.show(io::IO, nc::SampleClass)
  print(io, "----- Class ", nc.class, "-----\n")
  print(io, nc.members)
  return nothing
end
