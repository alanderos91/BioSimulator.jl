abstract type PriorityQueue end

struct PQBinaryHeap{K,V,O<:Ordering} <: PriorityQueue
  keys::Vector{K}
  values::Vector{V}
  keymap::Vector{Int} # map a key to its value
  valmap::Vector{Int} # map a value to its key

  function PQBinaryHeap{K,V,O}(o::O) where {K,V,O<:Ordering}
    new{K,V,O}(K[], V[], Int[], Int[])
  end

  function PQBinaryHeap{K,V,O}(keys, values) where {K,V,O<:Ordering}
    # sanity checks
    if length(keys) != length(values)
      throw(ArgumentError("number of keys and values is mismatched"))
    end

    if length(keys) != length(unique(keys))
      throw(ArgumentError("keys must be unique"))
    end

    # build association maps
    keymap = collect(1:length(keys))
    valmap = copy(keymap)
  
    # construct the object
    pq = new{K,V,O}(copy(keys), copy(values), keymap, valmap)

    # build the heap
    heapify!(pq)
  
    # obtain a sorting permutation for the keys array
    sortorder = sortperm(pq.keys)
  
    # do the sorting for the keys and keymap
    pq.keys   .= pq.keys[sortorder]
    pq.keymap .= pq.keymap[sortorder]

    # update the value to key associations
    for k in eachindex(sortorder)
      j = pq.keymap[k] # key k maps to position j in the heap,
      pq.valmap[j] = k # so position j must map back to k
    end
    
    return pq
  end
end

## heap traversal
@inline heap_left(j) = 2 * j
@inline heap_right(j) = 2 * j + 1
@inline heap_parent(j) = div(j, 2)

## low-level operations
@inline function sift_down!(pq::PQBinaryHeap, parent, bottom, is_violation)
  # unpack fields
  keys = pq.keys
  values = pq.values
  keymap = pq.keymap
  valmap = pq.valmap

  # save the parent value and start the search
  parent_value = values[parent]
  child = heap_left(parent)
  
  # while a left child node exists...
  while child <= bottom
    if child < bottom && is_violation(values[child], values[child + 1])
      child = child + 1
    end
    
    # if the parent does not violate the heap property, we're done
    if !is_violation(parent_value, values[child])
      break
    else
      node = heap_parent(child)
      
      node_keyix  = valmap[node]
      child_keyix = valmap[child]
      
      # perform swaps
      values[node] = values[child]
      
      keymap[node_keyix]  = child
      keymap[child_keyix] = node

      valmap[node]  = child_keyix
      valmap[child] = node_keyix
      
      child = heap_left(child)
    end
  end
  
  values[heap_parent(child)] = parent_value
  
  return pq
end

@inline function sift_up!(pq::PQBinaryHeap, top, child, is_violation)
  # unpack fields
  keys = pq.keys
  values = pq.values
  keymap = pq.keymap
  valmap = pq.valmap

  # while we're not at the top of the heap...
  while child > top
    parent = heap_parent(child)

    if is_violation(values[parent], values[child])
      parent_keyix = valmap[parent]
      child_keyix  = valmap[child]
      
      # swap values on the heap
      values[parent], values[child] = values[child], values[parent]
      
      keymap[parent_keyix] = child
      keymap[child_keyix]  = parent
      
      valmap[parent] = child_keyix
      valmap[child]  = parent_keyix
      
      # keep working up the tree
      child = parent
    else
      break
    end
  end
  
  return pq
end

function enqueue!(pq::PQBinaryHeap, value, key)
  # the value will be inserted at the bottom of the heap
  idx = length(pq) + 1
  
  # where should the key be inserted to maintain an ordered list?
  i = searchsortedfirst(pq.keys, key)

  insert!(pq.keys, i, key)   # insert `key` at index `i`
  insert!(pq.keymap, i, idx) # insert the association `i ---> idx`
  insert!(pq.valmap, idx, i) # insert the association `idx ---> i`
  push!(pq.values, value)    # add `value` at the end of the heap

  # fix upstream keys
  for j in i+1:idx
    pq.valmap[pq.keymap[j]] += 1
  end
  
  # check for heap order violations, working from the bottom to the top
  sift_up!(pq, 1, idx)
end

# TODO: make this efficient
function dequeue!(pq::PQBinaryHeap)
  # save the top value
  top_key, top_value = peektop(pq)

  # get associations for top and bottom nodes
  top_validx, bot_validx = 1, length(pq)
  top_keyidx, bot_keyidx = pq.valmap[top_validx], pq.valmap[bot_validx]

  # swap the nodes
  pq.values[top_validx], pq.values[bot_validx] = pq.values[bot_validx], pq.values[top_validx]

  pq.keymap[top_keyidx] = bot_validx
  pq.keymap[bot_keyidx] = top_validx
  pq.valmap[top_validx] = bot_keyidx
  pq.valmap[bot_validx] = top_keyidx

  if !isempty(pq)
    sift_down!(pq, 1, length(pq) - 1)
  end

  for i in length(pq):-1:top_keyidx+1
    pq.valmap[pq.keymap[i]] -= 1
  end

  pop!(pq.values)
  pop!(pq.valmap)
  deleteat!(pq.keys,   top_keyidx)
  deleteat!(pq.keymap, top_keyidx)

  # return the deleted entry
  return (top_key, top_value)
end

@inline function update_keyvalue!(pq::PQBinaryHeap, value, i)
  node = pq.keymap[i]
  old_value = pq.values[node]
  pq.values[node] = value

  return old_value
end

## heap property violation
is_minheap_violation(x, y) = y < x
is_maxheap_violation(x, y) = x < y

## type dispatch

### minheap
@inline function sift_down!(pq::PQBinaryHeap{K,V,O}, parent, bottom) where {K,V,O<:ForwardOrdering}
  sift_down!(pq, parent, bottom, is_minheap_violation)
end

@inline function sift_up!(pq::PQBinaryHeap{K,V,O}, top, child) where {K,V,O<:ForwardOrdering}
  sift_up!(pq, top, child, is_minheap_violation)
end

### maxheap
@inline function sift_down!(pq::PQBinaryHeap{K,V,O}, parent, bottom) where {K,V,O<:ReverseOrdering}
  sift_down!(pq, parent, bottom, is_maxheap_violation)
end

@inline function sift_up!(pq::PQBinaryHeap{K,V,O}, top, child) where {K,V,O<:ReverseOrdering}
  sift_up!(pq, top, child, is_maxheap_violation)
end

### heap construction
@inline function heapify!(pq::PQBinaryHeap)
  # a flat vector already enforces the shape property
  bottom = length(pq)

  # how many parent nodes are there?
  parents = heap_parent(bottom):-1:1

  # enforce the heap property at each parent node
  for parent in parents
    sift_down!(pq, parent, bottom)
  end
  
  return pq
end

## interface

# empty constructors
PQBinaryHeap{K,V}(o::Ordering=Forward) where {K,V} = PQBinaryHeap{K,V,typeof(o)}(o)

# non-empty constructor
PQBinaryHeap(keys, values, o::Ordering=Forward) = _pqbinheap_helper(o, keys, values)

_pqbinheap_helper(o::Ord, keys, values) where Ord = PQBinaryHeap{eltype(keys), eltype(values), Ord}(keys, values)

Base.length(pq::PQBinaryHeap)  = length(pq.keys)
Base.isempty(pq::PQBinaryHeap) = isempty(pq.keys)
Base.haskey(pq::PQBinaryHeap, key) = isempty(searchsorted(pq.keys, key)) ? false : true

peektop(pq::PQBinaryHeap) = (pq.keys[pq.valmap[1]], pq.values[1])

function Base.getindex(pq::PQBinaryHeap, key)
  irange = searchsorted(pq.keys, key)

  if isempty(irange)
    throw(KeyError(key))
  else
    return pq.values[pq.keymap[irange[end]]]
  end
end

function Base.get(pq::PQBinaryHeap, key, default)
  irange = searchsorted(pq.keys, key)

  return isempty(irange) ? default : pq.values[pq.keymap[irange[end]]]
end

function Base.setindex!(pq::PQBinaryHeap{K,V,ForwardOrdering}, value, key) where {K,V}
  _setindex!(pq, value, key, is_minheap_violation)
end

function Base.setindex!(pq::PQBinaryHeap{K,V,ReverseOrdering}, value, key) where {K,V}
  _setindex_kernel(pq, value, key, is_maxheap_violation)
end

@inline function _setindex!(pq, value, key, is_violation)
  i = get(pq, key, 0)

  # if the key exists...
  if i != 0
    old_value = update_keyvalue!(pq, value, i)
    
    if is_violation(value, old_value)
      sift_down!(pq, node, length(pq), is_violation)
    else
      sift_up!(pq, 1, node, is_violation)
    end
  # if they key does not exist...
  else
    enqueue!(pq, value, key)
  end
  
  return pq
end

function Base.iterate(pq::PQBinaryHeap, state = 1)
  state > length(pq) && return nothing

  key_idx = state
  val_idx = pq.keymap[key_idx]

  return ((pq.keys[key_idx], pq.values[val_idx]), state + 1)
end

Base.firstindex(pq::PQBinaryHeap) = 1
Base.lastindex(pq::PQBinaryHeap)  = length(pq)
