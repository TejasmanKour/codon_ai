def gc_energy(seq):
    """Approximate stability: GC stronger → lower ΔG"""
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")

    # simple energy model
    return -(2 * gc + 1 * at)


def sliding_delta_g(dna_seq, window=30):
    dna_seq = dna_seq.upper()
    values = []

    for i in range(len(dna_seq) - window + 1):
        sub = dna_seq[i:i+window]
        dg = gc_energy(sub)
        values.append(dg)

    return values


def detect_strong_structures(dg_values, threshold=-50):
    """
    Detect regions with strong secondary structure
    More negative ΔG = more stable = more problematic
    """
    strong_regions = []

    for i, val in enumerate(dg_values):
        if val < threshold:
            strong_regions.append({
                "position": i,
                "delta_g": val
            })

    return strong_regions