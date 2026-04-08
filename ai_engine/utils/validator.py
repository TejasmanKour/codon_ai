from Bio.Seq import Seq

def validate_sequence(seq):
    """
    Validate sequence using Biopython
    """
    try:
        Seq(seq)
        return True
    except Exception:
        return False


def detect_sequence_type(seq):
    """
    Detect whether input is DNA, Protein, or Invalid
    """
    if not seq:
        return "INVALID"

    seq = seq.upper().replace("\n", "").strip()

    # First validate using Biopython
    if not validate_sequence(seq):
        return "INVALID"

    dna_chars = set("ATGC")
    protein_chars = set("ACDEFGHIKLMNPQRSTVWY")

    seq_set = set(seq)

    # DNA check (strict)
    if seq_set.issubset(dna_chars):
        return "DNA"

    # Protein check
    elif seq_set.issubset(protein_chars):
        return "PROTEIN"

    # Mixed or invalid biological sequence
    else:
        return "INVALID"