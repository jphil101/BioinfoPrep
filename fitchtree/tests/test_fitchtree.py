import random
from itertools import combinations

import numpy as np
import pytest
from Bio.Phylo.TreeConstruction import DistanceMatrix

from fitchtree import fitch_margoliash, percent_standard_deviation

MOUNT = [[0, 22, 39, 39], [22, 0, 41, 41], [39, 41, 0, 18], [39, 41, 18, 0]]   # Mount, Bioinformatics


def leaf_lengths(tree):
    return {c.name: c.branch_length for c in tree.get_terminals()}


def path_lengths(tree, names):
    return {(x, y): tree.distance(x, y) for x, y in combinations(names, 2)}


def reference(rows, names):
    """Straightforward textbook version (recomputes every average from the taxa): leaf path lengths."""
    D = {(names[i], names[j]): rows[i][j] for i in range(len(names)) for j in range(len(names))}
    d = lambda X, Y: np.mean([D[i, j] for i in X for j in Y])
    C = {(n,): 0.0 for n in names}                      # cluster -> avg depth
    path = {n: {} for n in names}                       # leaf -> {cluster: distance from leaf up to cluster node}
    up = {(n,): {n: 0.0} for n in names}                # cluster -> {leaf: distance from cluster node down to leaf}
    while len(C) > 2:
        A, B = min(combinations(C, 2), key=lambda p: d(*p))
        ha, hb = C.pop(A), C.pop(B)
        R = sum(C, ())
        a = (d(A, B) + d(A, R) - d(B, R)) / 2
        b = d(A, B) - a
        up[A + B] = {**{l: v + a - ha for l, v in up[A].items()}, **{l: v + b - hb for l, v in up[B].items()}}
        C[A + B] = (len(A) * a + len(B) * b) / len(A + B)
    (X, hx), (Y, hy) = C.items()
    edge = d(X, Y) - hx - hy
    out = {}
    for x, y in combinations(names, 2):
        out[x, y] = path_between(x, y, up, edge, X, Y)
    return out


def path_between(x, y, up, edge, X, Y):
    # smallest cluster containing both, or the final edge
    common = [k for k in up if x in k and y in k]
    if not common:
        return up[X if x in X else Y][x] + edge + up[X if y in X else Y][y]
    k = min(common, key=len)
    child = [c for c in up if set(c) < set(k) and len(c) < len(k)]
    cx = max((c for c in child if x in c), key=len)
    cy = max((c for c in child if y in c), key=len)
    return (up[k][x] - up[cx][x]) + up[cx][x] + (up[k][y] - up[cy][y]) + up[cy][y]


def clock_like(n, seed):
    """Random ultrametric (clock-like, additive) distances on n taxa.
    Here the closest pair is always a true neighbour pair, so the method must recover the tree."""
    rng = random.Random(seed)
    groups = [([f"T{i}"], 0) for i in range(n)]          # (taxa, height)
    edges, height = [], 0
    while len(groups) > 1:
        (ga, ha), (gb, hb) = groups.pop(rng.randrange(len(groups))), groups.pop(rng.randrange(len(groups)))
        height += rng.randint(1, 9)
        edges += [(set(ga), height - ha), (set(gb), height - hb)]
        groups.append((ga + gb, height))
    names = [f"T{i}" for i in range(n)]
    rows = [[sum(k for e, k in edges if (x in e) != (y in e)) for y in names] for x in names]
    return rows, names


def test_mount_textbook_example():
    tree = fitch_margoliash(MOUNT, names="ABCD")
    assert leaf_lengths(tree) == {"A": 10, "B": 12, "C": 9, "D": 9}
    assert tree.distance("A", "C") == pytest.approx(10 + 20 + 9)          # internal branch = 20
    assert not tree.rooted


@pytest.mark.parametrize("n,seed", [(n, s) for n in range(3, 9) for s in range(5)])
def test_clock_like_distances_recovered_exactly(n, seed):
    rows, names = clock_like(n, seed)
    tree = fitch_margoliash(rows, names)
    for (x, y), v in path_lengths(tree, names).items():
        assert v == pytest.approx(rows[names.index(x)][names.index(y)])
    assert percent_standard_deviation(tree, rows, names) == pytest.approx(0, abs=1e-9)


@pytest.mark.parametrize("seed", range(30))
def test_matches_reference_on_random_non_additive(seed):
    rng = random.Random(seed)
    n = rng.randint(3, 7)
    names = [f"S{i}" for i in range(n)]
    rows = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        rows[i][j] = rows[j][i] = rng.randint(1, 20)
    tree = fitch_margoliash(rows, names)
    for pair, v in reference(rows, names).items():
        assert tree.distance(*pair) == pytest.approx(v)


def test_input_forms_agree():
    dm = DistanceMatrix(list("ABCD"), [r[:i + 1] for i, r in enumerate(MOUNT)])
    trees = [fitch_margoliash(dm), fitch_margoliash(MOUNT, "ABCD"), fitch_margoliash(np.array(MOUNT), list("ABCD"))]
    assert all(leaf_lengths(t) == leaf_lengths(trees[0]) for t in trees)


def test_default_names():
    assert sorted(leaf_lengths(fitch_margoliash(MOUNT))) == ["Seq1", "Seq2", "Seq3", "Seq4"]


def test_two_and_three_taxa():
    assert leaf_lengths(fitch_margoliash([[0, 6], [6, 0]], "AB")) == {"A": 3, "B": 3}
    assert leaf_lengths(fitch_margoliash([[0, 3, 4], [3, 0, 5], [4, 5, 0]], "ABC")) == {"A": 1, "B": 2, "C": 3}


def test_zero_distances_ok():
    tree = fitch_margoliash([[0, 0, 3, 5], [0, 0, 3, 5], [3, 3, 0, 2], [5, 5, 2, 0]], "ABCD")
    assert tree.distance("A", "B") == pytest.approx(0)
    assert percent_standard_deviation(tree, [[0, 0, 3, 5], [0, 0, 3, 5], [3, 3, 0, 2], [5, 5, 2, 0]], "ABCD") >= 0


def test_percent_standard_deviation_by_hand():
    # 3 taxa always fit exactly; check the formula on a fixed wrong tree instead
    tree = fitch_margoliash([[0, 2], [2, 0]], "AB")            # A-B path = 2
    assert percent_standard_deviation(tree, [[0, 4], [4, 0]], "AB") == pytest.approx(50)   # |4-2|/4 = 50%


def test_no_floating_point_noise():
    rows = [[0, 1, 3, 5, 6], [1, 0, 2, 4, 5], [3, 2, 0, 2, 3], [5, 4, 2, 0, 1], [6, 5, 3, 1, 0]]
    lengths = [c.branch_length for c in fitch_margoliash(rows).find_clades() if c.branch_length]
    assert all(abs(x) > 1e-9 for x in lengths)


def test_closest_pair_not_neighbours_is_a_known_limit():
    # additive tree ((T0:11,T2:1):5,(T1:7,T3:5)): T2,T3 are the closest pair (11) but not neighbours.
    # The clustering method joins them anyway (as UPGMA would); PSD reports the misfit.
    rows = [[0, 23, 12, 21], [23, 0, 13, 12], [12, 13, 0, 11], [21, 12, 11, 0]]
    tree = fitch_margoliash(rows)
    assert {c.name for c in tree.common_ancestor("Seq3", "Seq4").get_terminals()} == {"Seq3", "Seq4"}
    assert percent_standard_deviation(tree, rows) > 5


@pytest.mark.parametrize("rows,names,msg", [
    ([[0, 1], [1, 0, 2]], None, "square"),
    ([[0, 1], [2, 0]], None, "symmetric"),
    ([[0, 1], [1, 0]], ["A"], "unique names"),
    ([[0, 1], [1, 0]], ["A", "A"], "unique names"),
    ([[0, -1], [-1, 0]], None, "non-negative"),
    ([[0]], None, "at least 2"),
])
def test_bad_input(rows, names, msg):
    with pytest.raises(ValueError, match=msg):
        fitch_margoliash(rows, names)
