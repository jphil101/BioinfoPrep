try: import Bio, matplotlib, fitchtree
except ImportError: import os, sys; os.system(f'"{sys.executable}" -m pip install biopython matplotlib fitchtree')   # install if missing

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix
from fitchtree import fitch_margoliash

matrix = [[ 0, 22, 39, 39],
          [22,  0, 41, 41],
          [39, 41,  0, 18],
          [39, 41, 18,  0]]      # full matrix, or just the lower triangle
names = [f"Seq{i + 1}" for i in range(len(matrix))]
lower = []                                                         # DistanceMatrix wants the lower triangle
for i, row in enumerate(matrix):
    lower.append(row[:i + 1])
D = DistanceMatrix(names, lower)
tree = fitch_margoliash(D)


def branch_label(clade):                                           # text written on each branch
    if clade.branch_length:                                        # skip the root (None) and zero-length branches
        return round(clade.branch_length, 4)


Phylo.draw(tree, branch_labels=branch_label)
