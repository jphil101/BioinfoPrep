try: import Bio, matplotlib
except ImportError: import os, sys; os.system(f'"{sys.executable}" -m pip install biopython matplotlib')   # install if missing

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

n = int(input("How many sequences? "))
seqs = [input(f"Sequence {i + 1}: ").strip().upper() for i in range(n)]
names = [f"Seq{i + 1}" for i in range(n)]
mismatches = lambda s, t: sum(x != y for x, y in zip(s, t))        # distance = number of differing sites
D = DistanceMatrix(names, [[mismatches(s, t) for t in seqs[:i + 1]] for i, s in enumerate(seqs)])
print(D)                                                           # distance matrix
tree = DistanceTreeConstructor().nj(D)
Phylo.draw(tree, branch_labels=lambda c: c.branch_length and round(c.branch_length, 4))
