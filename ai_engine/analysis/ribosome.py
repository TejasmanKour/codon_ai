# Relative speed model (simplified)
FAST = 1.0
MEDIUM = 0.6
SLOW = 0.2

# Example classification
RARE_CODONS = {"AGG", "AGA", "ATA", "CTA", "CCC", "CGA"}

def codon_speed(codon):
    if codon in RARE_CODONS:
        return SLOW
    elif codon.endswith(("G", "C")):
        return FAST
    else:
        return MEDIUM


def ribosome_profile(dna_seq):
    dna_seq = dna_seq.upper()
    speeds = []

    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i:i+3]
        if len(codon) == 3:
            speeds.append(codon_speed(codon))

    return speeds