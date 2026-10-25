"""fitchtree -- Fitch-Margoliash distance trees for Biopython.

Builds a phylogenetic tree from a distance matrix with the Fitch-Margoliash
(1967) sequential clustering method, as taught in most bioinformatics courses:

1. join the closest pair of clusters A and B;
2. treat every other taxon as one composite group R and use the three-point
   formulas  a = (AB + AR - BR) / 2,  b = AB - a  for the branches to A and B;
3. subtract the average depth already inside A and B, merge, and repeat;
4. finish by joining the last three clusters at a single node (unrooted tree).

Distances between clusters are averages over their member taxa.

    >>> from fitchtree import fitch_margoliash
    >>> tree = fitch_margoliash([[0, 22, 39, 39],
    ...                          [22, 0, 41, 41],
    ...                          [39, 41, 0, 18],
    ...                          [39, 41, 18, 0]], names=["A", "B", "C", "D"])
    >>> print(tree.format("newick").strip())
    (A:10,B:12,(C:9,D:9):20):0;
"""

from Bio.Phylo.BaseTree import Clade, Tree
from Bio.Phylo.TreeConstruction import DistanceMatrix

__all__ = ["fitch_margoliash", "percent_standard_deviation"]
__version__ = "0.1.0"


def _to_distance_matrix(distances, names=None):
    """Accept a Biopython DistanceMatrix or a square matrix (list of lists / NumPy array)."""
    if isinstance(distances, DistanceMatrix):
        return distances
    rows = [[float(x) for x in row] for row in distances]
    n = len(rows)
    names = list(names) if names is not None else [f"Seq{i + 1}" for i in range(n)]
    if any(len(row) != n for row in rows):
        raise ValueError("distance matrix must be square")
    if len(names) != n or len(set(names)) != n:
        raise ValueError(f"need {n} unique names, got {names!r}")
    for i in range(n):
        for j in range(i):
            if abs(rows[i][j] - rows[j][i]) > 1e-9:
                raise ValueError(f"distance matrix is not symmetric at ({names[i]}, {names[j]})")
    return DistanceMatrix(names, [row[:i + 1] for i, row in enumerate(rows)])


def fitch_margoliash(distances, names=None):
    """Build an unrooted Fitch-Margoliash tree.

    Parameters
    ----------
    distances : Bio.Phylo.TreeConstruction.DistanceMatrix or square matrix
        Pairwise distances between taxa (list of lists or NumPy array also accepted).
    names : list of str, optional
        Taxon names when ``distances`` is a plain matrix (default Seq1, Seq2, ...).

    Returns
    -------
    Bio.Phylo.BaseTree.Tree
        Unrooted tree with branch lengths. Branch lengths can be negative when
        the distances are far from tree-like; that is a property of the method.
    """
    dm = _to_distance_matrix(distances, names)
    if len(dm.names) < 2:
        raise ValueError("need at least 2 taxa")
    if any(dm[i, j] < 0 for i in dm.names for j in dm.names):
        raise ValueError("distances must be non-negative")

    # each cluster: a key -> [Clade, number of taxa, average depth from its node to its taxa]
    clusters = {name: [Clade(name=name), 1, 0.0] for name in dm.names}
    dist = {a: {b: dm[a, b] for b in dm.names if b != a} for a in dm.names}   # average cluster distances

    def size(k):
        return clusters[k][1]

    def dist_to_rest(k, exclude):
        """Average distance from cluster k to every taxon outside k and `exclude`."""
        others = [x for x in clusters if x not in (k, exclude)]
        return sum(size(x) * dist[k][x] for x in others) / sum(size(x) for x in others)

    if len(clusters) == 2:
        a, b = dm.names
        for k in (a, b):
            clusters[k][0].branch_length = dm[a, b] / 2
        return _tree([clusters[a][0], clusters[b][0]])

    next_id = 0
    while len(clusters) > 3:
        keys = list(clusters)
        A, B = min(((x, y) for i, x in enumerate(keys) for y in keys[i + 1:]), key=lambda p: dist[p[0]][p[1]])
        dAB, dAR, dBR = dist[A][B], dist_to_rest(A, B), dist_to_rest(B, A)
        a = (dAB + dAR - dBR) / 2                      # node -> taxa of A (average)
        b = dAB - a                                    # node -> taxa of B (average)
        (ca, na, ha), (cb, nb, hb) = clusters.pop(A), clusters.pop(B)
        ca.branch_length, cb.branch_length = _clean(a - ha), _clean(b - hb)
        AB = f"_node{next_id}"
        next_id += 1
        dist[AB] = {}
        for x in clusters:                             # size-weighted average = average over taxa
            dist[AB][x] = dist[x][AB] = (na * dist[A][x] + nb * dist[B][x]) / (na + nb)
        for x in clusters:
            del dist[x][A], dist[x][B]
        del dist[A], dist[B]
        clusters[AB] = [Clade(clades=[ca, cb]), na + nb, (na * a + nb * b) / (na + nb)]

    A, B, C = clusters                                 # last three meet at one node
    for x, y, z in ((A, B, C), (B, A, C), (C, A, B)):
        clade, _, depth = clusters[x]
        clade.branch_length = _clean((dist[x][y] + dist[x][z] - dist[y][z]) / 2 - depth)
    return _tree([clusters[k][0] for k in (A, B, C)])


def percent_standard_deviation(tree, distances, names=None):
    """Fitch & Margoliash's goodness of fit (lower is better; 0 = tree fits exactly).

    PSD = 100 * sqrt( mean over pairs of ((observed - tree) / observed)^2 ).
    Pairs with an observed distance of 0 are skipped.
    """
    dm = _to_distance_matrix(distances, names)
    terms = [((dm[x, y] - tree.distance(x, y)) / dm[x, y]) ** 2
             for i, x in enumerate(dm.names) for y in dm.names[:i] if dm[x, y] != 0]
    return 100 * (sum(terms) / len(terms)) ** 0.5 if terms else 0.0


def _clean(x):
    """Drop floating-point noise like 2.8e-17."""
    return round(x, 12) + 0.0


def _tree(children):
    return Tree(root=Clade(clades=children), rooted=False)
