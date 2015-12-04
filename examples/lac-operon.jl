using BioSimulator
using DataFrames
using Gadfly

# String Constants
gI       = "gI"
rI       = "rI"
I        = "I"
o        = "o"
RNAp     = "RNAp"
r        = "r"
Z        = "Z"
lactose  = "lactose"
Ilactose = "Ilactose"
Io       = "Io"
RNApo    = "RNApo"

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
add_species!(network, gI,       1, false)
add_species!(network, rI,       0, false)
add_species!(network, I,       50)
add_species!(network, o,        1, false)
add_species!(network, RNAp,   100, false)
add_species!(network, r,        0)
add_species!(network, Z,        0)
add_species!(network, lactose, 20)
add_species!(network, Ilactose, 0, false)
add_species!(network, Io,       0, false)
add_species!(network, RNApo,    0, false)

# Reaction Definitions
add_reaction!(network, r1,  "c1");  set_parameter!(network, "c1",   c1)
add_reactant!(network, r1, gI)
add_product!( network, r1, gI)
add_product!( network, r1, rI)

add_reaction!(network, r2,  "c2");  set_parameter!(network, "c2",   c2)
add_reactant!(network, r2, rI)
add_product!( network, r2, rI)
add_product!( network, r2, I)

add_reaction!(network, r3,  "c3");  set_parameter!(network, "c3",   c3)
add_reactant!(network, r3, I)
add_reactant!(network, r3, lactose)
add_product!( network, r3, Ilactose)

add_reaction!(network, r4,  "c4");  set_parameter!(network, "c4",   c4)
add_reactant!(network, r4, Ilactose)
add_product!( network, r4, I)
add_product!( network, r4, lactose)

add_reaction!(network, r5,  "c5");  set_parameter!(network, "c5",   c5)
add_reactant!(network, r5, I)
add_reactant!(network, r5, o)
add_product!( network, r5, Io)

add_reaction!(network, r6,  "c6");  set_parameter!(network, "c6",   c6)
add_reactant!(network, r6, Io)
add_product!( network, r6, I)
add_product!( network, r6, o)

add_reaction!(network, r7,  "c7");  set_parameter!(network, "c7",   c7)
add_reactant!(network, r7, o)
add_reactant!(network, r7, RNAp)
add_product!( network, r7, RNApo)

add_reaction!(network, r8,  "c8");  set_parameter!(network, "c8",   c8)
add_reactant!(network, r8, RNApo)
add_product!( network, r8, o)
add_product!( network, r8, RNAp)

add_reaction!(network, r9,  "c9");  set_parameter!(network, "c9",   c9)
add_reactant!(network, r9, RNApo)
add_product!( network, r9, o)
add_product!( network, r9, RNAp)
add_product!( network, r9, r)

add_reaction!(network, r10, "c10"); set_parameter!(network, "c10", c10)
add_reactant!(network, r10, r)
add_product!( network, r10, r)
add_product!( network, r10, Z)

add_reaction!(network, r11, "c11"); set_parameter!(network, "c12", c11)
add_reactant!(network, r11, lactose)
add_reactant!(network, r11, Z)
add_product!( network, r11, Z)

add_reaction!(network, r12, "c12"); set_parameter!(network, "c12", c12)
add_reactant!(network, r12, rI)

add_reaction!(network, r13, "c13"); set_parameter!(network, "c13", c13)
add_reactant!(network, r13, I)

add_reaction!(network, r14, "c14"); set_parameter!(network, "c14", c14)
add_reactant!(network, r14, Ilactose)
add_product!( network, r14, lactose)

add_reaction!(network, r15, "c15"); set_parameter!(network, "c15", c15)
add_reactant!(network, r15, r)

add_reaction!(network, r16, "c16"); set_parameter!(network, "c16", c16)
add_reactant!(network, r16, Z)

@time df1 = simulate(network, tf=500.0, with=:ssa, dt=1.0, itr=1_000);
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
                         Species = repeat([lactose, r, Z, I], inner=[501]));

plot(new_df, x=:Time, y=:Mean, ymin=:Min, ymax=:Max, color=:Species, Geom.line)

df2  = df[ df[:Time] .== 500.0, :]
dfh = DataFrame()
for s in names(df2)
  if s == :Time; continue; end
  df = DataFrame(Time=df2[:Time], Count=df2[s], Species=s)
  dfh = vcat(dfh,df)
end

plot(dfh, x=:Count, Geom.density)
