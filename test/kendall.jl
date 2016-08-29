using BioSimulator
using Base.Test

function kendall_mean(i,t,α,μ,ν)
   x = exp((α-μ)*t)
  return i * x + ν/(α-μ)*(x-1)
end

X0 = 5
α  = 2.0
μ  = 1.0
ν  = 0.5

m = kendall(X0, α, μ, ν)

T = 4.0
dt = 0.1
seed = 5357
itr = 100_000

ssa = SSA(T)
frm = FRM(T)
nrm = NRM(T)
odm = ODM(T, 1000)
sal = SAL(T, 0.125, 100.0, 0.75)
algorithms = [ssa, frm, nrm, odm, sal]

# Compute mean for comparison
npts = round(Int, T / dt) + 1
t = linspace(0.0, T, npts)
theoretical = kendall_mean(X0,t,α,μ,ν)

# Run SSA and SAL once to compile
print("    Precompiling..."); @time begin
  for algorithm in algorithms
    simulate(m, algorithm, sampling_interval=dt, nrlz=1)
  end
end

print("    Running tests...\n\n")
for algorithm in algorithms
  print("   - Uniform ", uppercase(string(typeof(algorithm))))
  srand(seed); @time result = simulate(m, algorithm, sampling_interval=dt, nrlz=itr)
  computed = mean(result.data, 3)
  print("     |observed - theoretical| = ", abs(computed[end] - theoretical[end]), "\n")
  println()
end
