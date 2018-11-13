function F_autoreg(x, parms)
  gene, P2_gene, RNA, P, P2 = x
  k1, k1r, k2, k3, k4, k4r, k5, k6 = parms

  repression_binding = k1 * gene * P2
  reverse_repression = k1r * P2_gene
  transcription = k2 * gene
  translation = k3 * RNA
  dimerization = k4 * P * (P - 1) / 2
  dissociation = k4r * P2
  RNA_degradation = k5 * RNA
  protein_degradation = k6 * P

  return [repression_binding, reverse_repression, transcription, translation, dimerization, dissociation, RNA_degradation, protein_degradation]
end

function gillespie_autoreg(;
  gene = 10, P2_gene = 0, RNA = 0, P = 0, P2 = 0,
  k1 = 1.0, k1r = 10.0, k2 = 0.01, k3 = 10.0, k4 = 1.0, k4r = 1.0, k5 = 0.1, k6 = 0.01)
  # initial conditions
  x0 = [gene, P2_gene, RNA, P, P2]
  # stoichiometries
  nu = [
    [-1  1  0  0 -1]; # gene + P2 --> P2_gene
    [ 1 -1  0  0  1]; # P2_gene --> gene + P2
    [ 0  0  1  0  0]; # gene --> gene + RNA
    [ 0  0  0  1  0]; # RNA --> RNA + P
    [ 0  0  0 -2  1]; # P + P --> P2
    [ 0  0  0  2 -1]; # P2 --> P + P
    [ 0  0 -1  0  0]; # RNA --> 0
    [ 0  0  0 -1  0]  # P --> 0
  ]
  # parameter vector
  parms = [k1, k1r, k2, k3, k4, k4r, k5, k6]
  return x0, nu, parms
end

x0, nu, parms = gillespie_autoreg()
F = F_autoreg
