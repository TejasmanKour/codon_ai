import os
import json

# 🔥 Biopython import
from Bio.Seq import Seq

# 🔥 Import rare codon list
from ai_engine.analysis.rare_codon import RARE_CODONS


def load_table(org):
    base_dir = os.path.dirname(__file__)
    table_path = os.path.join(base_dir, "..", "codon_tables", f"{org}.json")
    table_path = os.path.abspath(table_path)

    with open(table_path) as f:
        return json.load(f)


# 🧠 Choose best codon (avoid rare codons)
def get_best_codon(codon_list):
    for codon in codon_list:
        if codon not in RARE_CODONS:
            return codon
    return codon_list[0]


# 🧬 NEW: Validate translation using Biopython
def validate_translation(dna_seq):
    try:
        protein = Seq(dna_seq).translate()
        return str(protein)
    except Exception:
        return None


def optimize_codon(sequence, organism):
    table = load_table(organism)
    dna = ""

    sequence = sequence.upper().strip()

    for aa in sequence:
        if aa in table:
            codons = table[aa]

            if isinstance(codons, list):
                best_codon = get_best_codon(codons)
            else:
                best_codon = codons

            dna += best_codon
        else:
            dna += "NNN"

    # 🧠 Validate final DNA translation
    translated = validate_translation(dna)

    if translated is None:
        # If something went wrong, return None (handled in app.py)
        return None

    return dna