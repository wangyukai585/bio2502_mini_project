# Mini Project: Dynamic Programming in Sequence Alignment

## 项目简介
本项目演示了两种经典的生物序列比对算法：
- **Needleman–Wunsch**（全局比对）
- **Smith–Waterman**（局部比对）

通过随机生成与人类 Chr1 相同碱基频率的 DNA 序列，并实施全局和局部比对，分析比对得分分布与示例对齐结果。

## 文件结构
```text
mini_project/
├── __pycache__/                    # Python 缓存目录
├── chr1.fa                         # 下载并解压后的 Chr1 序列文件
├── count_freqs.py                  # 计算碱基频率脚本
├── random_seqs.py                  # 按频率随机生成 DNA 序列脚本
├── compute_scores.py               # 全局比对得分矩阵计算脚本
├── plot_histograms.py              # 绘制直方图脚本
├── run_all.py                      # 串联全过程的主脚本
├── smith_waterman.py               # 面向对象的 Smith–Waterman 实现
├── local_alignment_demo.py         # 局部比对示例脚本（带回溯）
├── miniproject-dp.md               # 项目报告 Markdown 文件
├── requirements.txt                # Python 依赖列表
├── README.md                       # 本说明文件
├── hist_N50.png                    # N=50 全局比对得分直方图
├── hist_N100.png                   # N=100 全局比对得分直方图
├── hist_N200.png                   # N=200 全局比对得分直方图
├── hist_N500.png                   # N=500 全局比对得分直方图
├── scores_N50.npy                  # N=50 全局比对得分矩阵
├── scores_N100.npy                 # N=100 全局比对得分矩阵
├── scores_N200.npy                 # N=200 全局比对得分矩阵
├── scores_N500.npy                 # N=500 全局比对得分矩阵
├── local_scores_N50.npy            # N=50 局部比对得分矩阵示例