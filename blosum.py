try: import numpy, pandas
except ImportError: import os, sys; os.system(f'"{sys.executable}" -m pip install numpy pandas')   # install if missing

import numpy as np, pandas as pd
from itertools import combinations

n = int(input("How many sequences? "))
seqs = [input(f"Sequence {i + 1}: ").strip().upper() for i in range(n)]
x, y = zip(*[pair for col in zip(*seqs) for pair in combinations(col, 2)])   # pairs within each COLUMN
q = pd.crosstab(x + y, y + x, normalize=True)              # observed pair freq
f = q.sum(axis=1)                                          # residue freq
blosum = (2 * np.log2(q / np.outer(f, f))).round()         # 2 log2(observed / expected)
print(blosum)
