import pytest
from itertools import combinations
from Bio.Phylo.TreeConstruction import DistanceMatrix
from io import StringIO
from Bio import Phylo
from distances import distances
from upgma import upgma
from neighbor_joining import neighbor_joining
from fitch_margoliash import fitch_margoliash as fm_newick


def fitch_margoliash(D):
    return Phylo.read(StringIO(fm_newick(D)), "newick")


def matrix(names, rows):
    return DistanceMatrix(list(names), [r[:i + 1] for i, r in enumerate(rows)])


def leaf(tree, name):
    return next(tree.find_clades(name=name))


def siblings(tree, x, y):
    return {c.name for c in tree.common_ancestor(x, y).get_terminals()} == {x, y}


def assert_path_lengths_equal(tree, D):
    for x, y in combinations(D.names, 2):
        assert tree.distance(x, y) == pytest.approx(D[x, y])


# Wikipedia "UPGMA"/"WPGMA" worked example (5S rRNA)
WIKI_UPGMA = matrix("abcde", [[0, 17, 21, 31, 23], [17, 0, 30, 34, 21], [21, 30, 0, 28, 39],
                              [31, 34, 28, 0, 43], [23, 21, 39, 43, 0]])

# Wikipedia "Neighbor joining" worked example (additive)
WIKI_NJ = matrix("abcde", [[0, 5, 9, 9, 8], [5, 0, 10, 10, 9], [9, 10, 0, 8, 7],
                           [9, 10, 8, 0, 3], [8, 9, 7, 3, 0]])

# Mount, Bioinformatics, Fitch-Margoliash worked example
MOUNT_FM = matrix("ABCD", [[0, 22, 39, 39], [22, 0, 41, 41], [39, 41, 0, 18], [39, 41, 18, 0]])

SEQS = {"A": "ACGTACGTAC", "B": "ACGTACGTTC", "C": "ACGAACGATC", "D": "TCGAACCATC", "E": "TCGTTCCATG"}


def test_distances():
    D = distances(SEQS)
    assert D["A", "B"] == pytest.approx(0.1)   # 1 mismatch in 10
    assert D["A", "E"] == pytest.approx(0.6)


def test_upgma_is_wpgma_wikipedia():
    # Biopython's upgma averages the two merged clusters equally = WPGMA (the method taught in class)
    t = upgma(WIKI_UPGMA)
    assert siblings(t, "a", "b") and siblings(t, "c", "d")
    assert leaf(t, "a").branch_length == pytest.approx(8.5)
    assert leaf(t, "e").branch_length == pytest.approx(11)
    assert leaf(t, "c").branch_length == pytest.approx(14)
    assert all(t.distance(n) == pytest.approx(17.5) for n in "abcde")   # ultrametric, WPGMA root height 17.5


def test_nj_wikipedia():
    t = neighbor_joining(WIKI_NJ)
    assert {n: leaf(t, n).branch_length for n in "abcde"} == pytest.approx(dict(a=2, b=3, c=4, d=2, e=1))
    assert_path_lengths_equal(t, WIKI_NJ)                               # additive data recovered exactly


def test_fitch_mount():
    t = fitch_margoliash(MOUNT_FM)
    assert siblings(t, "A", "B") and siblings(t, "C", "D")
    assert {n: leaf(t, n).branch_length for n in "ABCD"} == pytest.approx(dict(A=10, B=12, C=9, D=9))
    assert t.distance("A", "C") - 10 - 9 == pytest.approx(20)           # internal branch = 20
    assert_path_lengths_equal(t, MOUNT_FM)


def test_fitch_additive_five_taxa():
    assert_path_lengths_equal(fitch_margoliash(WIKI_NJ), WIKI_NJ)


@pytest.mark.parametrize("method", [upgma, neighbor_joining, fitch_margoliash])
def test_from_sequences(method):
    t = method(distances(SEQS))
    assert sorted(c.name for c in t.get_terminals()) == list("ABCDE")
    assert siblings(t, "A", "B")                                         # closest pair grouped
