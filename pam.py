import numpy as np, pandas as pd
from itertools import combinations
np.seterr(divide="ignore")                       # log(0) = -inf, no warning


# PAM1 from closely related aligned seqs: count pairs between each PAIR OF SEQUENCES
def pam1(seqs):
    aa = sorted(set("".join(seqs)))
    F = pd.DataFrame(0.0, aa, aa)
    for s, t in combinations(seqs, 2):
        for x, y in zip(s, t):
            F.loc[x, y] += 1; F.loc[y, x] += 1
    f = F.sum(axis=1) / F.values.sum()             # residue freq
    M = F.div(F.sum(axis=1), axis=0)               # P(i -> j)
    I = np.eye(len(aa))
    return I + 0.01 / (f * (1 - np.diag(M))).sum() * (M - I), f   # scale to 1% change


def pam(seqs, n=250):
    M1, f = pam1(seqs)
    Mn = np.linalg.matrix_power(M1, n)             # PAMn = PAM1^n
    return pd.DataFrame(10 * np.log10(Mn / f.values), f.index, f.index).round()  # s = 10 log10(M_ij / f_j)


if __name__ == "__main__":
    print(pam(["ACDEAKLM", "ACDEAKLV", "ACEEAKLM", "SCDEAKIM"]))
