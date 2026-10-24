from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

seqs = {f"Seq{i + 1}": input(f"Sequence {i + 1}: ").strip().upper() for i in range(int(input("How many sequences? ")))}
names = list(seqs)
p = lambda s, t: sum(x != y for x, y in zip(s, t))                 # distance = number of mismatches
D = DistanceMatrix(names, [[p(seqs[a], seqs[b]) for b in names[:i + 1]] for i, a in enumerate(names)])
print(D)                                                           # distance matrix
tree = DistanceTreeConstructor().nj(D)
print(tree.format("newick"))                                       # tree as text, with branch lengths
Phylo.draw(tree, label_func=lambda c: c.name if c.is_terminal() else None,           # plot: leaf names only,
           branch_labels=lambda c: c.branch_length and f"{round(c.branch_length, 2):g}",   # branch lengths,
           yticks=([],), ylabel=("",))                                            # no y-axis clutter
