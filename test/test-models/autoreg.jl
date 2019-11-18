"""
- x0: [gene, P2_gene, RNA, P, P2]
- k1: repression binding rate
- k1r: reverse repression binding rate
- k2: transcription rate
- k3: translation rate
- k4: dimerization rate
- k4r: dissociation rate
- k5: RNA degradation
- k6: protein degradation rate
"""
function neg_autoreg(;x0 = [10, 0, 0, 0, 0], k1=1.0, k1r=10.0, k2=0.01, k3=10.0, k4=1.0, k4r=1.0, k5=0.1, k6=0.01)
  m = Network("auto-regulation")

  m <= Species("gene",    x0[1])
  m <= Species("P2_gene", x0[2])
  m <= Species("RNA",     x0[3])
  m <= Species("P",       x0[4])
  m <= Species("P2",      x0[5])

  m <= Reaction("repression binding",         k1,  "gene + P2 --> P2_gene")
  m <= Reaction("reverse repression binding", k1r, "P2_gene --> gene + P2")
  m <= Reaction("transcription",              k2,  "gene --> gene + RNA")
  m <= Reaction("translation",                k3,  "RNA --> RNA + P")
  m <= Reaction("dimerization",               k4,  "P + P --> P2")
  m <= Reaction("dissociation",               k4r, "P2 --> P + P")
  m <= Reaction("RNA degradation",            k5,  "RNA --> 0")
  m <= Reaction("protein degradation",        k6,  "P --> 0")

  return m
end
