# count_freqs.py
from collections import Counter

def count_bases(fasta_path):
    freqs = Counter()
    total = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq = line.strip().upper()
            freqs.update(seq)
            total += len(seq)
    # 归一化得到 A/C/G/T 的概率
    probs = {b: freqs[b]/total for b in "ACGT"}
    return probs

if __name__ == "__main__":
    p = count_bases("chr1.fa")
    print("A,C,G,T 频率:", p)
