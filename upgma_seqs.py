try: import Bio, matplotlib
except ImportError: import os, sys; os.system(f'"{sys.executable}" -m pip install biopython matplotlib')   # install if missing

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

seqs = ["ACGTACGTAC", "ACGTACGTTC", "ACGAACGATC", "TCGAACCATC", "TCGTTCCATG"]   # aligned, same length
names = [f"Seq{i + 1}" for i in range(len(seqs))]


def mismatches(s, t):                                              # distance = number of differing sites
    return sum(x != y for x, y in zip(s, t))


matrix = []                                                        # lower triangle of distances
for i, s in enumerate(seqs):
    row = []
    for t in seqs[:i + 1]:
        row.append(mismatches(s, t))
    matrix.append(row)
D = DistanceMatrix(names, matrix)
print(D)                                                           # distance matrix
tree = DistanceTreeConstructor().upgma(D)                          # (Biopython's upgma = WPGMA, as taught)


def branch_label(clade):                                           # text written on each branch
    if clade.branch_length:                                        # skip the root (None) and zero-length branches
        return round(clade.branch_length, 4)


Phylo.draw(tree, branch_labels=branch_label)
