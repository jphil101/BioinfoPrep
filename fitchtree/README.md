# fitchtree

Fitch–Margoliash distance trees for [Biopython](https://biopython.org).

Biopython builds UPGMA and neighbour-joining trees (`DistanceTreeConstructor`), but has no
Fitch–Margoliash method. `fitchtree` adds one, returning an ordinary `Bio.Phylo` tree that you can
draw, print or save like any other.

```bash
pip install fitchtree
```

## Usage

```python
from Bio import Phylo
from fitchtree import fitch_margoliash

tree = fitch_margoliash([[0, 22, 39, 39],
                         [22, 0, 41, 41],
                         [39, 41, 0, 18],
                         [39, 41, 18, 0]], names=["A", "B", "C", "D"])

print(tree.format("newick"))     # (A:10,B:12,(C:9,D:9):20):0;
Phylo.draw(tree, branch_labels=lambda c: c.branch_length and round(c.branch_length, 4))
```

`fitch_margoliash` accepts:

- a Biopython `DistanceMatrix` (e.g. from `DistanceCalculator`), or
- a square matrix (list of lists or NumPy array) plus optional `names` (default `Seq1`, `Seq2`, …).

`percent_standard_deviation(tree, distances, names=None)` gives Fitch & Margoliash's goodness-of-fit
score: 0 means the tree reproduces every distance exactly; lower is better.

## Method

This is the sequential clustering procedure from Fitch & Margoliash (1967), as presented in
textbooks such as Mount, *Bioinformatics: Sequence and Genome Analysis*:

1. Join the closest pair of clusters A and B.
2. Treat every other taxon as one group R and use the three-point formulas
   `a = (d(A,B) + d(A,R) − d(B,R)) / 2` and `b = d(A,B) − a` for the branches to A and B.
3. Subtract the average branch depth already inside A and B, merge them, and repeat.
4. Join the last three clusters at one central node (the tree is unrooted).

Distances between clusters are averages over their member taxa. When the distances fit a
clock-like tree, the method recovers it exactly (topology and branch lengths). Like UPGMA, it
always joins the closest pair first, so when the two closest taxa are not true neighbours (very
unequal rates of change) the topology can be wrong; `percent_standard_deviation` will show a
poor fit. Use neighbour joining for such data.

Branch lengths can come out negative when the distances are far from tree-like; that is a
property of the method, not an error.

This package does **not** perform the weighted least-squares search over topologies used by
PHYLIP's `fitch` program; it builds one tree by the clustering procedure above.

## License

MIT
