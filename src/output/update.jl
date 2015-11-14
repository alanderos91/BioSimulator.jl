function init_df(sname, tracked, n)
  x1 = [Float64]; x2 = [Int for i=1:length(tracked)]
  df = DataFrame(vcat(x1, x2), vcat([:Time], sname[tracked]), n)
end

function update!(::Explicit, df::DataFrame, n, t, t_next, dt, j, sname, tracked, spcs)
	insert_row!(df, t, tracked, spcs)
	return (t, j)
end

function update!(::Explicit, df::DataFrame, n, t_next, dt, j, sname, tracked, spcs)
	insert_row!(df, t_next, tracked, spcs)
	return (t_next, j)
end

function update!(::Uniform, df::DataFrame, n, t, t_next, dt, j, sname, tracked, spcs)
	while t >= t_next
		j = update_row!(df, t_next, sname, tracked, spcs, j)
		t_next = t_next + dt
		if j > n; return t_next, j; end
	end
	return (t_next, j)
end

function update!(::Uniform, df::DataFrame, n, t_next, dt, j, sname, tracked, spcs)
	while j <= n
		j = update_row!(df, t_next, sname, tracked, spcs, j)
		t_next = t_next + dt
	end
	return (t_next, j);
end

function insert_row!(df::DataFrame, t, tracked, spcs)
  x = spcs[tracked]
  push!(df, tuple(t, x...))
end

function update_row!(df::DataFrame, t, sname, tracked, spcs, j)
  df[j,1] = t
  for k in eachindex(tracked)
    df[j, 1+k] = spcs[tracked[k]]
  end
  return j + 1
end
