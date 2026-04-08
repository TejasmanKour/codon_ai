def gc_content_window(seq, window=10):
    seq = seq.upper()
    values = []

    for i in range(len(seq) - window + 1):
        sub = seq[i:i+window]
        gc = (sub.count("G") + sub.count("C")) / window
        values.append(gc)

    return values