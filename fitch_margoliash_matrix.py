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
D = DistanceMatrix(names, [row[:i + 1] for i, row in enumerate(matrix)])
tree = fitch_margoliash(D)
Phylo.draw(tree, branch_labels=lambda c: c.branch_length and round(c.branch_length, 4))
