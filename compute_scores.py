# compute_scores.py
import numpy as np
from Bio import pairwise2

def score_matrix(seqs, match=1, mismatch=-1, gap_open=-2, gap_extend=-0.5):
    M = len(seqs)
    mat = np.zeros((M, M))
    for i in range(M):
        for j in range(i, M):
            # globalms 使用 match/mismatch/gap 分数
            aln = pairwise2.align.globalms(
                seqs[i], seqs[j],
                match, mismatch,
                gap_open, gap_extend,
                score_only=True
            )
            mat[i,j] = mat[j,i] = aln
    return mat

if __name__ == "__main__":
    from random_seqs import generate_random_seqs
    from count_freqs import count_bases

    probs = count_bases("chr1.fa")
    for N in [50, 100, 200, 500]:
        seqs = generate_random_seqs(probs, N)
        print(f"正在计算 N={N} 的比对得分 …")
        mat = score_matrix(seqs)
        np.save(f"scores_N{N}.npy", mat)
