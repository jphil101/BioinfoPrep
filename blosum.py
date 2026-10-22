import numpy as np, pandas as pd
from itertools import combinations

def blosum(block):            # pairs within each COLUMN
    x, y = zip(*[p for col in zip(*block) for p in combinations(col, 2)])
    q = pd.crosstab(x + y, y + x, normalize=True)        # observed pair freq
    p = q.sum(axis=1)                                     # residue freq
    return (2 * np.log2(q / np.outer(p, p))).round()      # 2 log2(observed / expected)
