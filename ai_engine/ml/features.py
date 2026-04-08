import numpy as np

def extract_features(sequence):
    seq = sequence.upper()
    
    length = len(seq)
    gc = (seq.count('G') + seq.count('C')) / length if length > 0 else 0
    
    codon_count = len(seq) // 3
    
    return np.array([length, gc, codon_count]).reshape(1, -1)