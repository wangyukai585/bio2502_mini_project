# run_all.py
import os
from count_freqs import count_bases
from random_seqs import generate_random_seqs
from compute_scores import score_matrix
from plot_histograms import plot_hist
import numpy as np

def main():
    probs = count_bases("chr1.fa")
    for N in [50, 100, 200, 500]:
        print(f"\n=== N = {N} ===")
        seqs = generate_random_seqs(probs, N)
        mat = score_matrix(seqs)
        npy = f"scores_N{N}.npy"
        png = f"hist_N{N}.png"
        np.save(npy, mat)
        print(f"已保存 {npy}")
        plot_hist(mat, N)
        print(f"已保存 {png}")

if __name__ == "__main__":
    main()
