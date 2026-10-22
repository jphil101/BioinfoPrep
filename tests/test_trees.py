import pytest
from itertools import combinations
from conftest import tree_sequences


def leaf_lengths(tree):
    return {c.name: c.branch_length for c in tree.get_terminals()}


def siblings(tree, x, y):
    return {c.name for c in tree.common_ancestor(x, y).get_terminals()} == {x, y}


def assert_paths(tree, edges, L):
    """Every leaf-to-leaf path in the tree equals the true distance."""
    for x, y in combinations(sorted({n for e in edges for n in e}), 2):
        true = sum(k for e, k in edges.items() if (x in e) != (y in e))
        assert tree.distance(x, y) == pytest.approx(true / L)


def build(run, script, edges):
    lines, L = tree_sequences(edges)
    return run(script, str(len(lines)), *lines)["tree"], L


# ultrametric tree ((A:2,B:2):3,(C:4,D:4):1) -- UPGMA/WPGMA should recover it exactly
ULTRA = {"A": 2, "B": 2, "C": 4, "D": 4, "AB": 3, "CD": 1}
# Wikipedia "Neighbor joining" worked example: a2 b3 c4 d2 e1, internal 3 and 2
WIKI_NJ = {"a": 2, "b": 3, "c": 4, "d": 2, "e": 1, "ab": 3, "de": 2}
# Mount, Bioinformatics, Fitch-Margoliash worked example: A10 B12 C9 D9, internal 20
MOUNT_FM = {"A": 10, "B": 12, "C": 9, "D": 9, "AB": 20}


def test_distances_are_p_distances(run):
    g = run("upgma.py", "2", "A ACGTACGTAC", "B ACGTACGTTC")
    assert g["D"]["A", "B"] == pytest.approx(0.1)          # 1 mismatch in 10


def test_upgma_ultrametric(run):
    t, L = build(run, "upgma.py", ULTRA)
    assert siblings(t, "A", "B") and siblings(t, "C", "D")
    assert leaf_lengths(t) == pytest.approx({"A": 2 / L, "B": 2 / L, "C": 4 / L, "D": 4 / L})
    assert all(t.distance(n) == pytest.approx(5 / L) for n in "ABCD")   # all leaves at height 5
    assert_paths(t, ULTRA, L)


def test_upgma_is_wpgma(run):
    # Biopython's upgma averages the 2 merged clusters equally (WPGMA, as taught in class).
    # A,B close; C is 1 away from A,B... D far from C only: WPGMA vs UPGMA differ at the last merge.
    g = run("upgma.py", "4", "A AAAAAAAAAA", "B AAAAAAAAAT", "C AAAAAAATTT", "D TTTTTTTTTT")
    t = g["tree"]
    # merges: (A,B) h=.05 ; (AB,C) d=(.3+.2)/2=.25 -> h=.125 ; (ABC,D): WPGMA d=((1+.9)/2+.7)/2=.825 -> h=.4125
    assert t.distance("A") == pytest.approx(0.4125)      # true UPGMA would give (1+.9+.7)/3/2 = .4333


def test_nj_wikipedia(run):
    t, L = build(run, "neighbor_joining.py", WIKI_NJ)
    assert leaf_lengths(t) == pytest.approx({n: k / L for n, k in WIKI_NJ.items() if len(n) == 1})
    assert_paths(t, WIKI_NJ, L)


def test_fitch_mount(run):
    t, L = build(run, "fitch_margoliash.py", MOUNT_FM)
    assert siblings(t, "A", "B") and siblings(t, "C", "D")
    assert leaf_lengths(t) == pytest.approx({"A": 10 / L, "B": 12 / L, "C": 9 / L, "D": 9 / L})
    assert_paths(t, MOUNT_FM, L)                          # so internal branch = 20


def test_fitch_five_taxa(run):
    t, L = build(run, "fitch_margoliash.py", WIKI_NJ)
    assert_paths(t, WIKI_NJ, L)


@pytest.mark.parametrize("script", ["upgma.py", "neighbor_joining.py", "fitch_margoliash.py"])
def test_from_sequences(run, script):
    seqs = ["A ACGTACGTAC", "B ACGTACGTTC", "C ACGAACGATC", "D TCGAACCATC", "E TCGTTCCATG"]
    t = run(script, "5", *seqs)["tree"]
    assert sorted(c.name for c in t.get_terminals()) == list("ABCDE")
    assert siblings(t, "A", "B")                          # closest pair grouped
