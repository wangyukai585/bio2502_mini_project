# smith_waterman.py

import numpy as np

class SmithWaterman:
    """
    Smith–Waterman 本地比对实现。
    参数：
        match     : 匹配得分 (默认为 +2)
        mismatch  : 错配惩罚 (默认为 -1)
        gap       : 缺口惩罚 (默认为 -1)
    方法：
        score(seq1, seq2) -> (max_score, H, max_pos)
            计算得分矩阵 H 并返回最大得分及其位置。
        traceback(seq1, seq2, H, max_pos) -> (aligned1, aligned2)
            从 H 和 max_pos 回溯出一条最优局部比对。
    """
    def __init__(self, match=2, mismatch=-1, gap=-1):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap

    def score(self, seq1: str, seq2: str):
        m, n = len(seq1), len(seq2)
        H = np.zeros((m+1, n+1), dtype=int)
        max_score = 0
        max_pos = (0, 0)

        # 填表
        for i in range(1, m+1):
            for j in range(1, n+1):
                if seq1[i-1] == seq2[j-1]:
                    diag = H[i-1, j-1] + self.match
                else:
                    diag = H[i-1, j-1] + self.mismatch
                up   = H[i-1, j]   + self.gap
                left = H[i,   j-1] + self.gap
                H[i, j] = max(0, diag, up, left)
                if H[i, j] > max_score:
                    max_score = H[i, j]
                    max_pos   = (i, j)

        return max_score, H, max_pos

    def traceback(self, seq1: str, seq2: str, H: np.ndarray, max_pos):
        """
        从 max_pos 开始回溯，直到遇到 0 为止，重构一条最优局部比对串。
        返回 aligned1, aligned2（带 '-' 的对齐序列）。
        """
        aligned1, aligned2 = [], []
        i, j = max_pos
        while i > 0 and j > 0 and H[i, j] > 0:
            score_cur = H[i, j]
            if seq1[i-1] == seq2[j-1]:
                score_diag = H[i-1, j-1] + self.match
            else:
                score_diag = H[i-1, j-1] + self.mismatch

            if score_cur == score_diag:
                aligned1.append(seq1[i-1])
                aligned2.append(seq2[j-1])
                i -= 1; j -= 1
            elif score_cur == H[i-1, j] + self.gap:
                aligned1.append(seq1[i-1]); aligned2.append('-')
                i -= 1
            else:
                aligned1.append('-'); aligned2.append(seq2[j-1])
                j -= 1

        return ''.join(reversed(aligned1)), ''.join(reversed(aligned2))
