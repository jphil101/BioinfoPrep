from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

seqs = dict(input("Name and sequence (e.g. A ACGTAC): ").split() for _ in range(int(input("How many sequences? "))))
names = list(seqs)
p = lambda s, t: sum(x != y for x, y in zip(s, t)) / len(s)        # p-distance = fraction of mismatches
D = DistanceMatrix(names, [[p(seqs[a], seqs[b]) for b in names[:i + 1]] for i, a in enumerate(names)])
print(D)                                                           # distance matrix
tree = DistanceTreeConstructor().nj(D)
print(tree.format("newick"))                                       # tree as text, with branch lengths
Phylo.draw(tree, branch_labels=lambda c: c.branch_length and round(c.branch_length, 3))   # plot, branches labelled
