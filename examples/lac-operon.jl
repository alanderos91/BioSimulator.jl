using BioSimulator
using DataFrames
using Gadfly

r1  = "inhibitor transcription"
r2  = "inhibitor translation"
r3  = "lactose inhibitor binding"
r4  = "lactose inhibitor dissociation"
r5  = "inhibitor binding"
r6  = "inhibitor dissociation"
r7  = "RNAp binding"
r8  = "RNAp dissociation"
r9  = "transciption"
r10 = "translation"
r11 = "conversion"
r12 = "inhibitor RNA degradation"
r13 = "inhibitor degradation"
r14 = "lactose inhibitor degradation"
r15 = "RNA degradation"
r16 = "beta-galactosidase degradation"

# Parameter Values
c1  = 0.02
c2  = 0.1
c3  = 0.005
c4  = 0.1
c5  = 1.0
c6  = 0.01
c7  = 0.1
c8  = 0.01
c9  = 0.03
c10 = 0.1
c11 = 0.00005
c12 = 0.01
c13 = 0.002
c14 = 0.002
c15 = 0.01
c16 = 0.001

# Initialize network
network = Network("lac operon")

# Species Definitions
network <= Species(:gI,         1, istracked=false)
network <= Species(:rI,         0, istracked=false)
network <= Species(:I,         50)
network <= Species(:o,          1, istracked=false)
network <= Species(:RNAp,     100, istracked=false)
network <= Species(:r,          0)
network <= Species(:Z,          0)
network <= Species(:lactose,   20)
network <= Species(:Ilactose,   0, istracked=false)
network <= Species(:Io,         0, istracked=false)
network <= Species(:RNApo,      0, istracked=false)

# Reaction Definitions
network <= Reaction(r1,  :c1,  r=(:gI=>1),             p=(:gI=>1, :rI=>1))
network <= Reaction(r2,  :c2,  r=(:rI=>1),             p=(:rI=>1, :I=>1))
network <= Reaction(r3,  :c3,  r=(:I=>1, :lactose=>1), p=(:Ilactose=>1))
network <= Reaction(r4,  :c4,  r=(:Ilactose=>1),       p=(:lactose=>1, :I=>1))
network <= Reaction(r5,  :c5,  r=(:I=>1, :o=>1),       p=(:Io=>1))
network <= Reaction(r6,  :c6,  r=(:Io=>1),             p=(:I=>1, :o=>1))
network <= Reaction(r7,  :c7,  r=(:o=>1, :RNAp=>1),    p=(:RNApo=>1))
network <= Reaction(r8,  :c8,  r=(:RNApo=>1),          p=(:RNAp=>1, :o=>1))
network <= Reaction(r9,  :c9,  r=(:RNApo=>1),          p=(:RNAp=>1, :o=>1, :r=>1))
network <= Reaction(r10, :c10, r=(:r=>1),              p=(:r=>1, :Z=>1))
network <= Reaction(r11, :c11, r=(:lactose=>1, :Z=>1), p=(:Z=>1))
network <= Reaction(r12, :c12, r=(:rI=>1))
network <= Reaction(r13, :c13, r=(:I=>1))
network <= Reaction(r14, :c14, r=(:Ilactose=>1),       p=(:lactose=>1))
network <= Reaction(r15, :c15, r=(:r=>1))
network <= Reaction(r16, :c16, r=(:Z=>1))

network <= parameter(:c1,  c1)
network <= parameter(:c2,  c2)
network <= parameter(:c3,  c3)
network <= parameter(:c4,  c4)
network <= parameter(:c5,  c5)
network <= parameter(:c6,  c6)
network <= parameter(:c7,  c7)
network <= parameter(:c8,  c8)
network <= parameter(:c9,  c9)
network <= parameter(:c10, c10)
network <= parameter(:c11, c11)
network <= parameter(:c12, c12)
network <= parameter(:c13, c13)
network <= parameter(:c14, c14)
network <= parameter(:c15, c15)
network <= parameter(:c16, c16)

@time df1 = simulate(network, tf=500.0, with=:SSA, dt=1.0, itr=1_000);
      df = aggregate(df1, :Time, [mean,std]);

      new_df = DataFrame(Time    = repeat(convert(Vector, df[:Time]), outer=[4]),
                         Mean    = vcat(df[:lactose_mean], df[:r_mean], df[:Z_mean], df[:I_mean]),
                         Std     = vcat(df[:lactose_std],  df[:r_std],  df[:Z_std],  df[:I_std]),

                         Min     = vcat(df[:lactose_mean]  - df[:lactose_std],
                                       df[:r_mean]  - df[:r_std],
                                       df[:Z_mean] - df[:Z_std],
                                       df[:I_mean]  - df[:I_std]),

                         Max     = vcat(df[:lactose_mean] + df[:lactose_std],
                                       df[:r_mean]  + df[:r_std],
                                       df[:Z_mean] + df[:Z_std],
                                       df[:I_mean]  + df[:I_std]),
                         Species = repeat([:lactose, :r, :Z, :I], inner=[501]));

plot(new_df, x=:Time, y=:Mean, ymin=:Min, ymax=:Max, color=:Species, Geom.line)

# df2  = df[ df[:Time] .== 500.0, :]
# dfh = DataFrame()
# for s in names(df2)
#   if s == :Time; continue; end
#   df = DataFrame(Time=df2[:Time], Count=df2[s], Species=s)
#   dfh = vcat(dfh,df)
# end
#
# plot(dfh, x=:Count, Geom.density)
