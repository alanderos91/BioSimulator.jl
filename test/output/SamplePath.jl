import BioSimulator: build_output, update!

function make_sample_path(x_initial, tdata, states, save_points)
    xw = build_output(x_initial, nothing)

    for (t, state) in zip(tdata, states)
        update!(xw, t, state, save_points)
    end

    return xw
end

function check_sample_path(xw, ts, states)
    # always save initial state?
    @test length(xw) == length(ts) + 1

    # data should have been copied in the same order
    @assert length(states) == length(ts)

    for (i, state) in enumerate(states)
        @test xw.u[i+1] == state
        @test xw.t[i+1] == ts[i]
    end
end

@testset "SamplePath" begin
    n = 50
    d = 3
    tdata = sort!(10 * rand(n))
    states = [rand(0:n, d) for _ in tdata]

    # initialization
    x_initial = zeros(Int, d)

    xw = build_output(x_initial, nothing)
    # @test eltype(xw) == typeof(x_initial)

    @testset "save_points = nothing" begin
        xw = make_sample_path(x_initial, tdata, states, nothing)
        check_sample_path(xw, tdata, states)
    end

    @testset "save_points = iterable" begin
        stops_array = [0.5, 6.0]
        stops_range = 0.0:0.5:5
        # stops_tuple = (0.5, 6.0)

        # cases = (stops_array, stops_range, stops_tuple)
        cases = (stops_array, stops_range)

        @testset "$(Base.typename(typeof(stops)))" for stops in cases
            ts = filter(!isequal(0), stops)
            idxs = [searchsortedlast(tdata, t)+1 for t in ts]
            xw = make_sample_path(x_initial, tdata, states, stops)
            check_sample_path(xw, ts, states[idxs])
        end
    end
end
