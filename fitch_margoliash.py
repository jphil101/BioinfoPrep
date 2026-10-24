import numpy as np
from io import StringIO
from itertools import combinations
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix

n = int(input("How many sequences? "))
seqs = [input(f"Sequence {i + 1}: ").strip().upper() for i in range(n)]
names = [f"Seq{i + 1}" for i in range(n)]
mismatches = lambda s, t: sum(x != y for x, y in zip(s, t))        # distance = number of differing sites
D = DistanceMatrix(names, [[mismatches(s, t) for t in seqs[:i + 1]] for i, s in enumerate(seqs)])
print(D)                                                           # distance matrix

# join closest pair A,B; 3-point formula against R = everything else; repeat
d = lambda X, Y: np.mean([D[i, j] for i in X for j in Y])         # avg leaf-to-leaf distance
C = {(name,): (name, 0) for name in names}                         # cluster -> (newick, avg depth)
while len(C) > 2:
    A, B = min(combinations(C, 2), key=lambda pair: d(*pair))
    (na, ha), (nb, hb) = C.pop(A), C.pop(B)
    R = sum(C, ())                                                 # all other leaves
    a = (d(A, B) + d(A, R) - d(B, R)) / 2                          # branch to A
    b = d(A, B) - a                                                # branch to B
    C[A + B] = (f"({na}:{a - ha},{nb}:{b - hb})", (len(A) * a + len(B) * b) / len(A + B))
(X, (nx, hx)), (Y, (ny, hy)) = C.items()
tree = Phylo.read(StringIO(f"({nx},{ny}:{d(X, Y) - hx - hy});"), "newick")
Phylo.draw(tree, branch_labels=lambda c: c.branch_length and round(c.branch_length, 4))
