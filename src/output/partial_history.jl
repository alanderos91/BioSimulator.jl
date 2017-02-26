immutable SimData{T}
  data    :: T
  t_index :: LinSpace{Float64}
end

get_data(output :: SimData) = output.data

Base.getindex(output :: SimData, inds...) = getindex(output.data, inds...)
Base.setindex!(output :: SimData, inds...) = setindex!(output.data, inds...)

@inbounds function update!(
    output    :: SimData,
    Xt        :: Vector,
    t         :: Real,
    t_index   :: Range,
    epoch     :: Integer,
    trial     :: Integer
)
  data      = get_data(output)
  max_epoch = length(t_index)

  while (epoch <= max_epoch) && (t >= t_index[epoch])
    update!(data, Xt, epoch, trial)
    epoch = epoch + 1
  end

  return epoch
end

@inline function update!(
    output :: AbstractArray,
    Xt     :: Vector,
    epoch  :: Integer,
    trial  :: Integer
  )
  @inbounds for i in eachindex(Xt)
    output[i, epoch, trial] = Xt[i]
  end
end
