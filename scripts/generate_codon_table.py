import json
from collections import defaultdict
import os
import sys

# ==============================
# 🧬 Standard Genetic Code
# ==============================
CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",

    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",

    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",

    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}

# ==============================
# 📥 Read FASTA
# ==============================
def read_fasta(file_path):
    sequences = []
    with open(file_path) as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    sequences.append(seq)
                    seq = ""
            else:
                seq += line.strip().upper()
        if seq:
            sequences.append(seq)
    return sequences

# ==============================
# 🔢 Count Codons
# ==============================
def count_codons(sequences):
    codon_count = defaultdict(lambda: defaultdict(int))

    for seq in sequences:
        for i in range(0, len(seq)-2, 3):
            codon = seq[i:i+3]
            if len(codon) == 3 and codon in CODON_TABLE:
                aa = CODON_TABLE[codon]
                if aa != "*":  # skip stop codons
                    codon_count[aa][codon] += 1

    return codon_count

# ==============================
# 📊 Normalize (ranking)
# ==============================
def normalize_counts(codon_count):
    codon_usage = {}

    for aa in codon_count:
        sorted_codons = sorted(
            codon_count[aa].items(),
            key=lambda x: x[1],
            reverse=True
        )
        codon_usage[aa] = [codon for codon, count in sorted_codons]

    return codon_usage

# ==============================
# 🚀 MAIN FUNCTION
# ==============================
def main():
    base_dir = os.path.dirname(__file__)

    # ==============================
    # 📌 Get organism name from CLI
    # ==============================
    if len(sys.argv) < 2:
        print("❌ Usage: python generate_codon_table.py <organism>")
        print("Example: python generate_codon_table.py yeast")
        return

    organism = sys.argv[1].lower()

    # ==============================
    # 📂 Dynamic Paths
    # ==============================
    fasta_file = os.path.join(base_dir, "..", "data", f"{organism}_cds.fasta")
    output_file = os.path.join(base_dir, "..", "ai_engine", "codon_tables", f"{organism}.json")

    fasta_file = os.path.abspath(fasta_file)
    output_file = os.path.abspath(output_file)

    # ==============================
    # ❌ Check file exists
    # ==============================
    if not os.path.exists(fasta_file):
        print(f"❌ FASTA file not found: {fasta_file}")
        print("👉 Make sure file exists like: data/<organism>_cds.fasta")
        return

    # ==============================
    # 🔬 Processing
    # ==============================
    print(f"🧬 Processing organism: {organism}")
    print("Reading FASTA...")
    sequences = read_fasta(fasta_file)

    print("Counting codons...")
    codon_count = count_codons(sequences)

    print("Normalizing...")
    codon_usage = normalize_counts(codon_count)

    # ==============================
    # 💾 Save
    # ==============================
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, "w") as f:
        json.dump(codon_usage, f, indent=4)

    print(f"✅ Codon table saved to: {output_file}")

# ==============================
# ▶️ RUN
# ==============================
if __name__ == "__main__":
    main()