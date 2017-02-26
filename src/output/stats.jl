immutable MeanTrajectory{T}
  t_index  :: LinSpace{T}
  mean_mat :: Matrix{T}
  std_mat  :: Matrix{T}
end

function MeanTrajectory(output :: SimData)
  data    = get_data(output)
  t_index = get_time(output)

  n, m, _  = size(data)

  mean_mat = transpose(reshape(mean(data, 3), m, n))
  std_mat  = transpose(reshape(std(data, 3), m, n))

  return MeanTrajectory(t_index, mean_mat, std_mat)
end
