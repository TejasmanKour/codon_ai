from collections import Counter

def codon_usage(dna_seq):
    dna_seq = dna_seq.upper()
    codons = [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3) if len(dna_seq[i:i+3]) == 3]

    return Counter(codons)