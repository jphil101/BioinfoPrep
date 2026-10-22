import numpy as np, pandas as pd
from itertools import combinations

block = [input(f"Sequence {i + 1}: ") for i in range(int(input("How many sequences? ")))]
x, y = zip(*[pair for col in zip(*block) for pair in combinations(col, 2)])   # pairs within each COLUMN
q = pd.crosstab(x + y, y + x, normalize=True)              # observed pair freq
p = q.sum(axis=1)                                          # residue freq
blosum = (2 * np.log2(q / np.outer(p, p))).round()         # 2 log2(observed / expected)
print(blosum)
