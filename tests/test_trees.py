import pytest
import matplotlib.pyplot as plt
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


def named(tree, leaves):
    """Scripts name sequences Seq1, Seq2, ...; rename to the test's leaf letters."""
    for c in tree.get_terminals():
        c.name = leaves[int(c.name[3:]) - 1]
    return tree


def build(run, script, edges):
    seqs, leaves = tree_sequences(edges)
    return named(run(script, str(len(seqs)), *seqs)["tree"], leaves), 1   # distances are mismatch counts


# ultrametric tree ((A:2,B:2):3,(C:4,D:4):1) -- UPGMA/WPGMA should recover it exactly
ULTRA = {"A": 2, "B": 2, "C": 4, "D": 4, "AB": 3, "CD": 1}
# Wikipedia "Neighbor joining" worked example: a2 b3 c4 d2 e1, internal 3 and 2
WIKI_NJ = {"a": 2, "b": 3, "c": 4, "d": 2, "e": 1, "ab": 3, "de": 2}
# Mount, Bioinformatics, Fitch-Margoliash worked example: A10 B12 C9 D9, internal 20
MOUNT_FM = {"A": 10, "B": 12, "C": 9, "D": 9, "AB": 20}


def test_distances_are_mismatch_counts(run):
    g = run("upgma.py", "2", "ACGTACGTAC", "  acgtacgttc ")    # stray spaces / lowercase are fine
    assert g["D"]["Seq1", "Seq2"] == 1                            # 1 mismatch


def test_upgma_ultrametric(run):
    t, L = build(run, "upgma.py", ULTRA)
    assert siblings(t, "A", "B") and siblings(t, "C", "D")
    assert leaf_lengths(t) == pytest.approx({"A": 2 / L, "B": 2 / L, "C": 4 / L, "D": 4 / L})
    assert all(t.distance(n) == pytest.approx(5 / L) for n in "ABCD")   # all leaves at height 5
    assert_paths(t, ULTRA, L)


def test_upgma_is_wpgma(run):
    # Biopython's upgma averages the 2 merged clusters equally (WPGMA, as taught in class).
    # A,B close; C is 1 away from A,B... D far from C only: WPGMA vs UPGMA differ at the last merge.
    g = run("upgma.py", "4", "AAAAAAAAAA", "AAAAAAAAAT", "AAAAAAATTT", "TTTTTTTTTT")
    t = named(g["tree"], "ABCD")
    # merges: (A,B) h=.5 ; (AB,C) d=(3+2)/2=2.5 -> h=1.25 ; (ABC,D): WPGMA d=((10+9)/2+7)/2=8.25 -> h=4.125
    assert t.distance("A") == pytest.approx(4.125)       # true UPGMA would give (10+9+7)/3/2 = 4.333


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
    seqs = ["ACGTACGTAC", "ACGTACGTTC", "ACGAACGATC", "TCGAACCATC", "TCGTTCCATG"]
    t = named(run(script, "5", *seqs)["tree"], "ABCDE")
    assert sorted(c.name for c in t.get_terminals()) == list("ABCDE")
    assert siblings(t, "A", "B")                          # closest pair grouped


@pytest.mark.parametrize("script", ["upgma.py", "neighbor_joining.py", "fitch_margoliash.py"])
def test_output_shows_distances_and_branch_lengths(run, script):
    plt.close("all")
    g = run(script, "4", "ACGTACGTAC", "ACGTACGTTC", "ACGAACGATC", "TCGAACCATC")
    tree, printed = g["tree"], g["_printed"]
    assert str(g["D"]) in printed                                             # distance matrix printed
    labels = {t.get_text() for t in plt.gca().texts}                          # plot: every branch labelled
    for c in tree.find_clades():
        if round(c.branch_length or 0, 4):
            assert str(round(c.branch_length, 4)) in labels
    assert {"Seq1", "Seq2", "Seq3", "Seq4"} <= {l.strip() for l in labels}     # leaves labelled
