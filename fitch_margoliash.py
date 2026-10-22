import numpy as np
from itertools import combinations

def fitch_margoliash(D):      # join closest pair, 3-point formula vs. the rest, repeat
    d = lambda X, Y: np.mean([D[i, j] for i in X for j in Y])     # avg leaf-to-leaf distance
    C = {(n,): (n, 0) for n in D.names}                            # cluster -> (newick, avg depth)
    while len(C) > 2:
        A, B = min(combinations(C, 2), key=lambda p: d(*p))
        (na, ha), (nb, hb) = C.pop(A), C.pop(B)
        R = sum(C, ())                                             # all other leaves
        a = (d(A, B) + d(A, R) - d(B, R)) / 2
        b = d(A, B) - a
        C[A + B] = (f"({na}:{a-ha},{nb}:{b-hb})", (len(A)*a + len(B)*b) / len(A+B))
    (X, (nx, hx)), (Y, (ny, hy)) = C.items()
    return f"({nx},{ny}:{d(X, Y)-hx-hy});"                         # Newick tree
