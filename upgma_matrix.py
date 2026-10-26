try: import Bio, matplotlib
except ImportError: import os, sys; os.system(f'"{sys.executable}" -m pip install biopython matplotlib')   # install if missing

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

matrix = [[ 0, 17, 21, 31, 23],
          [17,  0, 30, 34, 21],
          [21, 30,  0, 28, 39],
          [31, 34, 28,  0, 43],
          [23, 21, 39, 43,  0]]      # full matrix, or just the lower triangle
names = [f"Seq{i + 1}" for i in range(len(matrix))]
D = DistanceMatrix(names, [row[:i + 1] for i, row in enumerate(matrix)])
tree = DistanceTreeConstructor().upgma(D)                          # (Biopython's upgma = WPGMA, as taught)
Phylo.draw(tree, branch_labels=lambda c: c.branch_length and round(c.branch_length, 4))
