# Reactions
RepressionBinding:
    P2 + Gene > P2Gene
    RepressionBinding_k1*Gene*P2

ProteinDegradation:
    P > $pool
    ProteinDegradation_k6 * P

Dimerisation:
    {2.0}P > P2
    Dimerisation_k4 * 0.5 * P * (P-1)

Dissociation:
    P2 > {2} P
    Dissociation_k4r * P2

Transcription:
    $pool > Rna
    Transcription_k2 * Gene

RnaDegradation:
    Rna > $pool
    RnaDegradation_k5 * Rna

Translation:
    $pool > P
    Translation_k3 * Rna

ReverseRepressionBinding:
    P2Gene > P2 + Gene
    ReverseRepressionBinding_k1r * P2Gene
 
# Fixed species
 
# Variable species
P2 = 0.0
P = 0.0
Rna = 0.0
Gene = 10.0
P2Gene = 0.0
 
# Parameters
RepressionBinding_k1 = 1.0
ProteinDegradation_k6 = 0.01
Dimerisation_k4 = 1.0
Dissociation_k4r = 1.0
Transcription_k2 = 0.01
RnaDegradation_k5 = 0.1
Translation_k3 = 10.0
ReverseRepressionBinding_k1r = 10.0 
