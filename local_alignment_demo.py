# local_alignment_demo.py

import numpy as np
from count_freqs    import count_bases
from random_seqs    import generate_random_seqs
from smith_waterman import SmithWaterman

def compute_local_matrix(seqs, aligner):
    M = len(seqs)
    mat = np.zeros((M, M), dtype=int)
    for i in range(M):
        for j in range(i, M):
            score, _, _ = aligner.score(seqs[i], seqs[j])
            mat[i, j] = mat[j, i] = score
    return mat

if __name__ == "__main__":
    # 示例：N=50, 序列数 M=10（为了快速演示）
    N, M = 50, 10
    probs = count_bases("chr1.fa")
    seqs  = generate_random_seqs(probs, N, M)
    aligner = SmithWaterman(match=2, mismatch=-1, gap=-1)

    # 1) 计算局部比对得分矩阵
    local_mat = compute_local_matrix(seqs, aligner)
    print("Local score matrix shape:", local_mat.shape)
    print(local_mat)

    # 2) 对第 0 条和第 1 条序列回溯一次
    s1, s2 = seqs[0], seqs[1]
    max_score, H, max_pos = aligner.score(s1, s2)
    aln1, aln2 = aligner.traceback(s1, s2, H, max_pos)
    print("\nSeq[0]:", s1)
    print("Seq[1]:", s2)
    print("Local alignment score:", max_score)
    print("Alignment result:")
    print(aln1)
    print(aln2)

    # 3) 保存矩阵到文件
    np.save(f"local_scores_N{N}.npy", local_mat)
    print(f"已保存 local_scores_N{N}.npy")
