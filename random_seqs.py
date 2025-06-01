# random_seqs.py
import random

def generate_random_seqs(probs, N, M=1000):
    """返回 M 条长度为 N、基于给定频率随机生成的序列列表"""
    bases, weights = zip(*probs.items())
    seqs = [
        "".join(random.choices(bases, weights, k=N))
        for _ in range(M)
    ]
    return seqs

if __name__ == "__main__":
    from count_freqs import count_bases
    probs = count_bases("chr1.fa")
    # 示例：N=50 的随机序列
    seqs50 = generate_random_seqs(probs, 50)
