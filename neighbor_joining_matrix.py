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
D = DistanceMatrix(names, [row[:i + 1] for i, row in enumerate(matrix)])
tree = DistanceTreeConstructor().nj(D)
Phylo.draw(tree, branch_labels=lambda c: c.branch_length and round(c.branch_length, 4))
