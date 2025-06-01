# plot_histograms.py
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Heiti TC'

def plot_hist(mat, N):
    idx = np.triu_indices_from(mat, k=1)
    scores = mat[idx]
    plt.figure()
    plt.hist(scores, bins=50)
    plt.title(f"全局比对得分分布 (N={N})")
    plt.xlabel("得分")
    plt.ylabel("频数")
    plt.grid(True)
    plt.savefig(f"hist_N{N}.png")
    plt.close()

if __name__ == "__main__":
    for N in [50, 100, 200, 500]:
        mat = np.load(f"scores_N{N}.npy")
        plot_hist(mat, N)
        print(f"已保存 hist_N{N}.png")
