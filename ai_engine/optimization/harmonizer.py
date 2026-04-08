import os
import json
import random


def load_table(org):
    base_dir = os.path.dirname(__file__)
    table_path = os.path.join(base_dir, "..", "codon_tables", f"{org}.json")
    table_path = os.path.abspath(table_path)

    with open(table_path) as f:
        return json.load(f)


def classify_codons(codon_list):
    """
    Splits codons into high / medium / low frequency groups
    Assumes codons are ordered by frequency
    """
    n = len(codon_list)

    if n == 1:
        return {"high": codon_list, "medium": [], "low": []}

    return {
        "high": codon_list[: max(1, n // 3)],
        "medium": codon_list[max(1, n // 3): max(2, 2 * n // 3)],
        "low": codon_list[max(2, 2 * n // 3):],
    }


def harmonize_codon(sequence, organism):
    table = load_table(organism)
    dna = ""

    sequence = sequence.upper().strip()

    for aa in sequence:
        if aa in table:
            codons = table[aa]

            if isinstance(codons, list):
                groups = classify_codons(codons)

                # 🎯 Randomly choose based on natural distribution
                choice_type = random.choices(
                    ["high", "medium", "low"],
                    weights=[0.6, 0.3, 0.1]
                )[0]

                selected_group = groups[choice_type]

                if selected_group:
                    dna += random.choice(selected_group)
                else:
                    dna += codons[0]
            else:
                dna += codons
        else:
            dna += "NNN"

    return dna