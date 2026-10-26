try: import Bio, matplotlib
except ImportError: import os, sys; os.system(f'"{sys.executable}" -m pip install biopython matplotlib')   # install if missing

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

matrix = [[ 0,  5,  9,  9,  8],
          [ 5,  0, 10, 10,  9],
          [ 9, 10,  0,  8,  7],
          [ 9, 10,  8,  0,  3],
          [ 8,  9,  7,  3,  0]]      # full matrix, or just the lower triangle
names = [f"Seq{i + 1}" for i in range(len(matrix))]
lower = []                                                         # DistanceMatrix wants the lower triangle
for i, row in enumerate(matrix):
    lower.append(row[:i + 1])
D = DistanceMatrix(names, lower)
tree = DistanceTreeConstructor().nj(D)


def branch_label(clade):                                           # text written on each branch
    if clade.branch_length:                                        # skip the root (None) and zero-length branches
        return round(clade.branch_length, 4)


Phylo.draw(tree, branch_labels=branch_label)
