import BioSimulator: build_output, update_samplepath!, save_state, get_regular_path

function make_sample_path(f, x_initial, tdata, states, save_points)
    xw = build_output(f, nothing, x_initial, nothing, nothing)

    for (t, state) in zip(tdata, states)
        update_samplepath!(f, xw, nothing, t, state, nothing, save_points)
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

function save_ones(simulator, state, model)
    return ones(length(state))
end

@testset "SamplePath" begin
    n = 50
    d = 3
    tdata = sort!(10 * rand(n))
    states = [rand(0:n, d) for _ in tdata]

    # initialization
    x_initial = zeros(Int, d)

    @testset "save_points = nothing" begin
        xw = make_sample_path(save_state, x_initial, tdata, states, nothing)
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
            xw = make_sample_path(save_state, x_initial, tdata, states, stops)
            check_sample_path(xw, ts, states[idxs])
        end
    end

    @testset "save_function" begin
        xw = make_sample_path(save_ones, x_initial, tdata, states, nothing)
        # not happy with this test
        expected = save_ones.(nothing, states, nothing)
        check_sample_path(xw, tdata, expected)
    end

    @testset "regular paths" begin
        stops = 0.0:0.5:5
        ts = filter(!isequal(0), stops)
        idxs = [searchsortedlast(tdata, t)+1 for t in ts]

        # make SamplePath from full observations
        xw = make_sample_path(save_state, x_initial, tdata, states, nothing)

        # make SamplePath from partial observations
        yw = make_sample_path(save_state, x_initial, tdata, states, stops)

        # transform full SamplePath to 'regular' path
        zw = get_regular_path(xw, stops)

        @test yw == zw
    end

    @testset "simulate" begin
        network = Network("extinction test")
        network <= Species("X", 1)
        network <= Reaction("birth", 1e-6, "X --> X + X")
        network <= Reaction("death", 1e2, "X --> 0")

        # case: save_points = nothing
        # expected behavior: tfinal is included
        trajectory = simulate(network, Direct(), tfinal = 10.0)

        @test trajectory.t[end] == 10.0
    end
end
