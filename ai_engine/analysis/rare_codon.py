# Example: E. coli rare codons (can expand per organism)
RARE_CODONS = {
    "AGG", "AGA", "ATA", "CTA", "CCC", "CGA"
}

def find_rare_codons(dna_seq):
    dna_seq = dna_seq.upper()
    rare_positions = []

    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i:i+3]
        if codon in RARE_CODONS:
            rare_positions.append((i, codon))

    return rare_positions