using BioSimulator

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
network <= Species("gI",         1)
network <= Species("rI",         0)
network <= Species("I",         50)
network <= Species("o",          1)
network <= Species("RNAp",     100)
network <= Species("r",          0)
network <= Species("Z",          0)
network <= Species("lactose",   20)
network <= Species("Ilactose",   0)
network <= Species("Io",         0)
network <= Species("RNApo",      0)

# Reaction Definitions
network <= Reaction(r1,  c1,  :(gI --> gI + rI))
network <= Reaction(r2,  c2,  :(rI --> rI + I))
network <= Reaction(r3,  c3,  :(I + lactose --> Ilactose))
network <= Reaction(r4,  c4,  :(Ilactose --> I + lactose))
network <= Reaction(r5,  c5,  :(I + o --> Io))
network <= Reaction(r6,  c6,  :(Io --> I + o))
network <= Reaction(r7,  c7,  :(o + RNAp --> RNApo))
network <= Reaction(r8,  c8,  :(RNApo --> o + RNAp))
network <= Reaction(r9,  c9,  :(RNApo --> o + RNAp + r))
network <= Reaction(r10, c10, :(r --> r + Z))
network <= Reaction(r11, c11, :(lactose + Z --> Z)
network <= Reaction(r12, c12, :(rI --> 0))
network <= Reaction(r13, c13, :(I --> 0))
network <= Reaction(r14, c14, :(Ilactose --> lactose))
network <= Reaction(r15, c15, :(r --> 0))
network <= Reaction(r16, c16, :(Z --> 0))

@time result = simulate(network, SSA(500.0), sampling_interval=1.0, nrlz=1_000)

plot(MeanTrajectory([:lactose, :r, :Z, :I]))
