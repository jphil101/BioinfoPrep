import numpy as np, pandas as pd
from itertools import combinations

n = int(input("How many sequences? "))
seqs = [input(f"Sequence {i + 1}: ").strip().upper() for i in range(n)]
pam_number = int(input("PAM number (e.g. 250): "))
x, y = zip(*[pair for s, t in combinations(seqs, 2) for pair in zip(s, t)])   # pairs between each PAIR OF SEQUENCES
q = pd.crosstab(x + y, y + x, normalize=True)              # observed pair freq
f = q.sum(axis=1)                                          # residue freq
M = q.div(f, axis=0)                                       # P(i -> j) = q_ij / f_i
I = np.eye(len(f))
M1 = I + 0.01 / (1 - np.trace(q)) * (M - I)                # scale so 1% of residues change = PAM1
Mn = pd.DataFrame(np.linalg.matrix_power(M1, pam_number), f.index, f.index)   # PAMn = PAM1^n
pam = (10 * np.log10(Mn / f)).round()                      # 10 log10(M_ij / f_j)
print(pam)
