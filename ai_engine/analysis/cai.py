def calculate_cai(dna_seq):
    dna_seq = dna_seq.upper()
    codons = [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3)]

    score = 0
    valid_codons = 0

    for codon in codons:
        if len(codon) == 3:
            # Simple scoring (can improve later)
            if codon.endswith(("G", "C")):
                score += 1
            valid_codons += 1

    if valid_codons == 0:
        return 0

    return round(score / valid_codons, 3)