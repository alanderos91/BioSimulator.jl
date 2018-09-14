using DataFrames

x0 = zeros(Int, 5)

x0[1] = 9 # gene = 9
x0[2] = 1 # P2_gene = 1
x0[3] = 2 # RNA = 2
x0[4] = 5 # P = 5
x0[5] = 3 # P2 = 3

model = autoreg(x0 = x0)

srand(5357)

# single simulation
result = simulate(model, Direct(), Val(:full), time = 100.0)
mydf = DataFrame(result)
mydf_columns = names(mydf)

# make sure we have the correct column names
@test :trial ∈ mydf_columns
@test :time ∈ mydf_columns
@test :P ∈ mydf_columns
@test :gene ∈ mydf_columns
@test :RNA ∈ mydf_columns
@test :P2_gene ∈ mydf_columns
@test :P2 ∈ mydf_columns

# check trial number and time columns
@test all(n -> n == 1, mydf[:trial])
@test mydf[:time][1] == 0.0
@test all(t -> t < 100.0, mydf[:time])

# check initial values for each population
@test mydf[:gene][1]    == x0[1]
@test mydf[:P2_gene][1] == x0[2]
@test mydf[:RNA][1]     == x0[3]
@test mydf[:P][1]       == x0[4]
@test mydf[:P2][1]      == x0[5]

# check properties of each population
@test all(x -> x == 10, mydf[:gene] + mydf[:P2_gene]) # 10 copies
@test all(x -> x < 10, mydf[:RNA]) # RNA population should be small
