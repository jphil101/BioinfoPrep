import numpy as np, pandas as pd
from itertools import combinations

def pam(seqs, n=250):         # pairs between each PAIR OF SEQUENCES
    x, y = zip(*[p for s, t in combinations(seqs, 2) for p in zip(s, t)])
    q = pd.crosstab(x + y, y + x, normalize=True)        # observed pair freq
    f = q.sum(axis=1)                                     # residue freq
    M = q.div(f, axis=0) - np.eye(len(f))                 # P(i -> j) - I
    M1 = np.eye(len(f)) + 0.01 / (1 - np.trace(q)) * M   # scale to 1% change = PAM1
    Mn = pd.DataFrame(np.linalg.matrix_power(M1, n), f.index, f.index)   # PAMn = PAM1^n
    return (10 * np.log10(Mn / f)).round()                # 10 log10(M_ij / f_j)
