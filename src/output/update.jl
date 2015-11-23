function init_df(sname, tracked, n)
  x1 = [Float64]; x2 = [Int for i=1:length(tracked)]
  df = DataFrame(vcat(x1, x2), vcat([:Time], sname[tracked]), n)
  for s in sname
    df[s] = zeros(Int, size(df,1))
  end
  return df
end

#immutable Output
#  df::DataFrame
#  meta::Dict{ASCIIString,Any}
#end

type Updater
  dt::Float64
  t_next::Float64

  tracked::Vector{Int}

  blocksize::Int
  maxindex::Int
  index::Int
end

function Updater(dt, tracked, blocksize)
  return Updater(dt, 0.0, tracked, blocksize, blocksize, 1)
end

#update!(otype,  o, u, species)
#update!(otype, df, u, species)

# Safeguard against unsafe operation?
function update!(::Explicit, df, t, u, species)
  x = species[u.tracked]
  push!(df, tuple(t, x...))
  u.index += 1
end

function update!(::Uniform, df, t, u, species)
  tracked = u.tracked
  while t >= u.t_next
    if u.index > u.maxindex; return; end
    i = u.index

    df[i,1] = u.t_next
    for j in eachindex(tracked)
      df[i,1+j] = species[tracked[j]]
    end

    u.t_next += u.dt
    u.index += 1
  end
end

function update!(::Mean, df, t, u, species)
  tracked = u.tracked
  while t >= u.t_next
    i = u.index;

    df[i,1] = u.t_next
    for j in eachindex(tracked)
      df[i,1+j] += species[tracked[j]]
    end

    u.t_next += u.dt
    u.index += 1
  end
end

function update!(::Histogram, df, t, u, species)
  return;
end

function final_update!(::Explicit, df, t, u, species)
  x = species[u.tracked]
  push!(df, tuple(t, x...))
  u.index += 1
end

function final_update!(::Uniform, df, t, u, species)
  tracked = u.tracked
  while u.index <= u.maxindex
    i = u.index;

    df[i,1] = u.t_next
    for j in eachindex(tracked)
      df[i,1+j] = species[tracked[j]]
    end

    u.t_next += u.dt
    u.index += 1
  end

  u.t_next = 0.0
  if u.maxindex < size(df,1)
    u.maxindex = u.maxindex + u.blocksize
  end
end

function final_update!(::Mean, df, t, u, species)
  tracked = u.tracked
  while u.index <= u.maxindex
    i = u.index;

    df[i,1] = u.t_next
    for j in eachindex(tracked)
      df[i,1+j] += species[tracked[j]]
    end

    u.t_next += u.dt
    u.index += 1
  end

  u.t_next = 0.0
  u.index = 1
end

function final_update!(::Histogram, df, t, u, species)
  i = u.index; tracked = u.tracked

  df[i,1] = t
  for j in eachindex(tracked)
    df[i,1+j] = species[tracked[j]]
  end

  u.index += 1
end
