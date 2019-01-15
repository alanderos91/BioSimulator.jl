"""
  Prokaryotic auto-regulation network
"""
function autoreg(;
  gene = 10, P2_gene = 0, RNA = 0, P = 0, P2 = 0,
  k1  :: Float64 = 0.01,
  k2  :: Float64 = 10.0,
  k3  :: Float64 = 1.0,
  k3r :: Float64 = 1.0,
  k4  :: Float64 = 1.0,
  k4r :: Float64 = 10.0,
  k5  :: Float64 = 0.1,
  k6  :: Float64 = 0.01)
  m = Network("auto-regulation")

  m <= Species("gene",    gene)
  m <= Species("P2_gene", P2_gene)
  m <= Species("RNA",     RNA)
  m <= Species("P",       P)
  m <= Species("P2",      P2)

  m <= Reaction("transcription",              k1,  "gene --> gene + RNA")
  m <= Reaction("translation",                k2,  "RNA --> RNA + P")
  m <= Reaction("dimerization",               k3,  "P + P --> P2")
  m <= Reaction("dissociation",               k3r, "P2 --> P + P")
  m <= Reaction("repression binding",         k4,  "gene + P2 --> P2_gene")
  m <= Reaction("reverse repression binding", k4r, "P2_gene --> gene + P2")
  m <= Reaction("RNA degradation",            k5,  "RNA --> 0")
  m <= Reaction("protein degradation",        k6,  "P --> 0")

  return m
end

model = autoreg()
