import numpy as np, pandas as pd
from itertools import combinations
np.seterr(divide="ignore")                       # log(0) = -inf, no warning


# BLOSUM from an ungapped block: count pairs within each COLUMN
def blosum(block):
    aa = sorted(set("".join(block)))
    F = pd.DataFrame(0.0, aa, aa)
    for col in zip(*block):
        for x, y in combinations(col, 2):
            F.loc[x, y] += 1; F.loc[y, x] += 1
    q = F / F.values.sum()                         # observed pair freq
    p = q.sum(axis=1)                              # residue freq
    return (2 * np.log2(q / np.outer(p, p))).round()   # s = 2 log2(observed / expected)


if __name__ == "__main__":
    print(blosum(["ACDA", "ACEA", "SCDA", "ACDG"]))
