from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

seqs = dict(input("Name and sequence (e.g. A ACGTAC): ").split() for _ in range(int(input("How many sequences? "))))
names = list(seqs)
p = lambda s, t: sum(x != y for x, y in zip(s, t)) / len(s)        # p-distance = fraction of mismatches
D = DistanceMatrix(names, [[p(seqs[a], seqs[b]) for b in names[:i + 1]] for i, a in enumerate(names)])
tree = DistanceTreeConstructor().nj(D)
Phylo.draw_ascii(tree)
