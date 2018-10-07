using Random

import Base.Order: ForwardOrdering, ReverseOrdering, Forward, Reverse
import NewBioSimulator: enqueue!, dequeue!, heapify!

@testset "PQBinaryHeap" begin
  # Test dequeing in sorted order.
  function test_issorted!(pq::PQBinaryHeap, rev=false)
    (_, last_val) = dequeue!(pq)
    while !isempty(pq)
      (_, next_val) = dequeue!(pq)
      if !rev
        @test last_val <= next_val
      else
        @test next_val <= last_val
      end
      last_val = next_val
    end
  end

  # TODO...
  function test_isrequested!(pq::PQBinaryHeap, keys)
    i = 0
    while !isempty(pq)
      krqst    = keys[i+=1]
      krcvd, _ = dequeue!(pq, krqst)
      @test krcvd == krqst
    end
  end

  # create the PQBinaryHeap but do not build the heap
  # function init_maxpq(keys, values)
  #   K = eltype(keys)
  #   V = eltype(values)

  #   keymap = collect(1:length(keys))
  #   valmap = copy(keymap)

  #   return PQBinaryHeap{K,V,ForwardOrdering}(keys, values, keymap, valmap)
  # end

  # function init_minpq(keys, values)
  #   K = eltype(keys)
  #   V = eltype(values)
    
  #   keymap = collect(1:length(keys))
  #   valmap = copy(keymap)

  #   return PQBinaryHeap{K,V,ReverseOrdering}(keys, values, keymap, valmap)
  # end

  # create the PQBinaryHeap by building the heap
  function build_maxpq(keys, values)
    K = eltype(keys)
    V = eltype(values)

    return PQBinaryHeap{K,V,ReverseOrdering}(keys, values)
  end

  function build_minpq(keys, values)
    K = eltype(keys)
    V = eltype(values)

    return PQBinaryHeap{K,V,ForwardOrdering}(keys, values)
  end

  @testset "Low-Level Heap Operations" begin
    # pmax = 1000
    # n = 10000

    # keys, values = collect(1:n), rand(1:pmax, n)

    # @testset "heapify!" begin
    #   maxpq = init_maxpq(keys, values)
    #   minpq = init_minpq(keys, values)
      
    #   heapify!(maxpq)
    #   heapify!(minpq)

    #   @test issorted([dequeue!(maxpq)[2] for _ in length(values)])
    #   @test issorted([dequeue!(minpq)[2] for _ in length(values)])
    # end

    # TODO: What else should be in here?
  end

  @testset "Constructors" begin
    pmax = 1000
    n = 10000

    keys, values = collect(1:n), rand(1:pmax, n)

    @testset "build from parallel arrays" begin
      test_issorted!(build_minpq(keys, values), false)
      test_issorted!(build_maxpq(keys, values), true)
    end

    @testset "different array sizes => ArgumentError" begin
      @test_throws ArgumentError PQBinaryHeap(['a'], [1, 2])
      @test_throws ArgumentError PQBinaryHeap(['a', 'b'], [1])
    end

    @testset "duplicate key => ArgumentError" begin
      @test_throws ArgumentError PQBinaryHeap(['a', 'a'], [1, 2])
    end
  end

  @testset "Methods" begin
    pmax = 1000
    n = 10000

    easy_sorted_keys,   easy_sorted_values   = ['a', 'b'],   [1, 2]
    easy_unsorted_keys, easy_unsorted_values = ['b', 'a'],   [2, 1]
    
    easy_cases = [
      "sorted/sorted"     => (easy_sorted_keys,   easy_sorted_values),
      "unsorted/sorted"   => (easy_unsorted_keys, easy_sorted_values),
      "sorted/unsorted"   => (easy_sorted_keys,   easy_unsorted_values),
      "unsorted/unsorted" => (easy_unsorted_keys, easy_unsorted_values)
    ]

    random_unsorted_keys, random_unsorted_values = shuffle(1:n),               rand(1:pmax, n)
    random_sorted_keys,   random_sorted_values   = sort(random_unsorted_keys), sort(random_unsorted_values)

    random_cases = [
      "sorted/sorted"     => (random_sorted_keys,   random_sorted_values),
      "unsorted/sorted"   => (random_unsorted_keys, random_sorted_values),
      "sorted/unsorted"   => (random_sorted_keys,   random_unsorted_values),
      "unsorted/unsorted" => (random_unsorted_keys, random_unsorted_values)
    ]

    examples = ["easy" => easy_cases, "random" => random_cases]
    # examples = ["easy" => easy_cases]

    @testset "$(exampletype) cases" for (exampletype, example) in examples
      @testset "$(case)" for (case, problem) in example
        (keys, values) = problem

        minpq = build_minpq(keys, values)
        maxpq = build_maxpq(keys, values)

        # we can only check for the correct value, key mappings will be stored in a non-deterministic fashion
        @testset "peektop" begin
          minval, minidx = findmin(values)
          maxval, maxidx = findmax(values)

          @test peektop(minpq)[2] == minval
          @test peektop(maxpq)[2] == maxval
        end

        @testset "get/getindex" begin
          if exampletype == "easy"
            # key does not exist
            @test get(minpq, 'c', 0) == 0
            @test get(maxpq, 'c', 0) == 0

            @test_throws KeyError minpq['c']
            @test_throws KeyError maxpq['c']

            # key exists
            a_value = values[something(findfirst(isequal('a'), keys), 0)]

            @test get(minpq, 'a', 0) == a_value
            @test get(maxpq, 'a', 0) == a_value
            @test minpq['a'] == a_value
            @test maxpq['a'] == a_value
          else
            # key does not exist
            @test get(minpq, n + 1, 0) == 0
            @test get(maxpq, n + 1, 0) == 0

            @test_throws KeyError minpq[n + 1]
            @test_throws KeyError maxpq[n + 1]

            # key exists
            random_key   = rand(keys)
            random_value = values[something(findfirst(isequal(random_key), keys), 0)]

            @test get(minpq, random_key, 0) == random_value
            @test get(maxpq, random_key, 0) == random_value
            @test minpq[random_key] == random_value
            @test maxpq[random_key] == random_value
          end
        end

        @testset "enqueue!" begin
          K, V = eltype(keys), eltype(values)

          minpq_copy = PQBinaryHeap{K,V}(Forward)
          maxpq_copy = PQBinaryHeap{K,V}(Reverse)

          for i in eachindex(keys)
            key   = keys[i]
            value = values[i]

            enqueue!(minpq_copy, value, key)
            enqueue!(maxpq_copy, value, key)
          end

          test_issorted!(minpq_copy, false)
          test_issorted!(maxpq_copy, true)
        end

        # again, we can only test the value returned by dequeue!
        @testset "dequeue!" begin
          minpq_copy = deepcopy(minpq)
          maxpq_copy = deepcopy(maxpq)
          
          sorted_order = sortperm(values)

          for i in eachindex(sorted_order)
            ifwd = sorted_order[i]
            irev = sorted_order[end - i + 1]

            nextmin = values[ifwd]
            nextmax = values[irev]

            @test dequeue!(minpq_copy)[2] == nextmin
            @test dequeue!(maxpq_copy)[2] == nextmax
          end
        end

        @testset "setindex!" begin
          @testset "changing priorities" begin
            if exampletype == "easy"
              minpq_copy = deepcopy(minpq)
              maxpq_copy = deepcopy(maxpq)

              minpq_copy['a'] = 100
              maxpq_copy['a'] = 100

              # make sure we do not add a duplicate key
              @test length(minpq_copy) == length(minpq)
              @test length(maxpq_copy) == length(maxpq)

              # was the value updated correctly?
              @test minpq_copy['a'] == 100
              @test maxpq_copy['a'] == 100
            end
          end

          @testset "adding a new item" begin
            
          end
        end

        # remember, key mappings are non-deterministic
        @testset "Iteration" begin
          K, V = eltype(keys), eltype(values)

          minpq_copy = PQBinaryHeap{K,V}(Forward)
          maxpq_copy = PQBinaryHeap{K,V}(Reverse)

          for (key, value) in minpq
            enqueue!(minpq_copy, value, key)
          end

          @test minpq.keys   == minpq_copy.keys
          @test minpq.values[minpq.keymap] == minpq_copy.values[minpq_copy.keymap]

          for (key, value) in maxpq
            enqueue!(maxpq_copy, value, key)
          end

          @test maxpq.keys == maxpq_copy.keys
          @test maxpq.values[maxpq.keymap] == maxpq_copy.values[maxpq_copy.keymap]
        end
      end
    end
  end
end
