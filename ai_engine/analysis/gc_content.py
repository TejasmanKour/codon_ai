from Bio.SeqUtils import gc_fraction

def gc_content(seq):
    """
    Calculate GC content using Biopython (returns percentage)
    """
    if not seq:
        return 0

    seq = seq.upper().strip()

    try:
        return gc_fraction(seq) * 100
    except Exception:
        return 0