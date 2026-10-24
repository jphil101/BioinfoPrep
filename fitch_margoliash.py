import numpy as np
from io import StringIO
from itertools import combinations
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix

seqs = {f"Seq{i + 1}": input(f"Sequence {i + 1}: ").strip().upper() for i in range(int(input("How many sequences? ")))}
names = list(seqs)
p = lambda s, t: sum(x != y for x, y in zip(s, t))                 # distance = number of mismatches
D = DistanceMatrix(names, [[p(seqs[a], seqs[b]) for b in names[:i + 1]] for i, a in enumerate(names)])
print(D)                                                           # distance matrix

# join closest pair, 3-point formula vs. the rest, repeat
d = lambda X, Y: np.mean([D[i, j] for i in X for j in Y])         # avg leaf-to-leaf distance
C = {(n,): (n, 0) for n in names}                                  # cluster -> (newick, avg depth)
while len(C) > 2:
    A, B = min(combinations(C, 2), key=lambda pair: d(*pair))
    (na, ha), (nb, hb) = C.pop(A), C.pop(B)
    R = sum(C, ())                                                 # all other leaves
    a = (d(A, B) + d(A, R) - d(B, R)) / 2
    b = d(A, B) - a
    C[A + B] = (f"({na}:{a - ha},{nb}:{b - hb})", (len(A) * a + len(B) * b) / len(A + B))
(X, (nx, hx)), (Y, (ny, hy)) = C.items()
tree = Phylo.read(StringIO(f"({nx},{ny}:{d(X, Y) - hx - hy});"), "newick")
print(tree.format("newick"))                                       # tree as text, with branch lengths
Phylo.draw(tree, label_func=lambda c: c.name if c.is_terminal() else None,           # plot: leaf names only,
           branch_labels=lambda c: c.branch_length and f"{round(c.branch_length, 2):g}",   # branch lengths,
           yticks=([],), ylabel=("",))                                            # no y-axis clutter
