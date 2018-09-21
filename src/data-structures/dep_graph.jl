# make sure we have the required imports
import Base: iterate, IteratorSize, IteratorEltype, eltype, length

abstract type DependencyGraph end

@inline dependents(dep_graph::DependencyGraph, j) = dependents(dg_iterator(dep_graph), dep_graph, j)

abstract type DGIterStyle end

struct DGView <: DGIterStyle end
struct DGLazy <: DGIterStyle end

# Implementation: range tuples + single partitioned vector

struct DGVector{Iter<:DGIterStyle} <: DependencyGraph
  rxn_range::Vector{Tuple{Int,Int}}
  rxn_deps::Vector{Int}
end

## constructor 

function DGVector{Iter}(model::Network) where Iter <: DGIterStyle
  num_reactions = number_reactions(model)
  reactions = reaction_list(model)
  dg = model.dep_graph

  indexmap = OrderedDict(key => i for (i, key) in enumerate(keys(reactions)))

  rxn_range = Vector{Tuple{Int,Int}}(undef, num_reactions)
  rxn_deps  = Int[]

  idx = 1
  for (key, j) in indexmap
    start = idx
    for dep in dg[key]
      push!(rxn_deps, indexmap[dep])
      idx += 1
    end
    stop = idx - 1

    rxn_range[j] = (start, stop)
  end

  return DGVector{Iter}(rxn_range, rxn_deps)
end

## iteration using views

@inline function dependents_view(dg::DGVector, j)
  @inbounds (start, stop) = dg.rxn_range[j]
  
  return view(dg.rxn_deps, start:stop)
end

## iteration using laziness

struct DGVectorIterator
  deps::Vector{Int}
  start::Int
  stop::Int
end

Base.IteratorSize(::DGVectorIterator) = Base.HasLength()
Base.IteratorEltype(::DGVectorIterator) = Base.HasEltype()

Base.length(iter::DGVectorIterator) = iter.stop - iter.start + 1
Base.eltype(iter::DGVectorIterator) = eltype(iter.deps)

Base.iterate(iter::DGVectorIterator) = iterate(iter, iter.start)

function Base.iterate(iter::DGVectorIterator, state)
  return state > iter.stop ? nothing : (iter.deps[state], state + 1)
end

@inline function dependents_lazy(dg::DGVector, j)
  @inbounds (start, stop) = dg.rxn_range[j]
  
  return DGVectorIterator(dg.rxn_deps, start, stop)
end

## iteration protocol

@inline dg_iterator(dg::DGVector{Iter}) where Iter <: DGIterStyle = Iter()

@inline dependents(::DGView, dep_graph::DGVector, j) = dependents_view(dep_graph, j)
@inline dependents(::DGLazy, dep_graph::DGVector, j) = dependents_lazy(dep_graph, j)
