import pytest
import matplotlib.pyplot as plt
from itertools import combinations
from conftest import tree_sequences

METHODS = ["upgma", "neighbor_joining", "fitch_margoliash"]
SCRIPTS = [f"{m}_{kind}.py" for m in METHODS for kind in ("seqs", "matrix")]


def leaf_lengths(tree):
    return {c.name: c.branch_length for c in tree.get_terminals()}


def siblings(tree, x, y):
    """x, y are neighbours: some branch separates exactly {x, y} from the rest (works rooted or unrooted)."""
    leaves = {c.name for c in tree.get_terminals()}
    splits = [{c.name for c in clade.get_terminals()} for clade in tree.find_clades()]
    return any(s == {x, y} or leaves - s == {x, y} for s in splits)


def true_distance(edges, x, y):
    return sum(k for e, k in edges.items() if (x in e) != (y in e))


def assert_paths(tree, edges):
    """Every leaf-to-leaf path in the tree equals the true distance."""
    for x, y in combinations(sorted({n for e in edges for n in e}), 2):
        assert tree.distance(x, y) == pytest.approx(true_distance(edges, x, y))


def named(tree, leaves):
    """Scripts name taxa Seq1, Seq2, ...; rename to the test's leaf letters."""
    for c in tree.get_terminals():
        c.name = leaves[int(c.name[3:]) - 1]
    return tree


def from_seqs(run_with, method, edges):
    seqs, leaves = tree_sequences(edges)
    return named(run_with(f"{method}_seqs.py", seqs=seqs)["tree"], leaves)


def from_matrix(run_with, method, edges):
    leaves = sorted({n for e in edges for n in e})
    matrix = [[true_distance(edges, x, y) for y in leaves] for x in leaves]
    return named(run_with(f"{method}_matrix.py", matrix=matrix)["tree"], leaves)


# ultrametric tree ((A:2,B:2):3,(C:4,D:4):1) -- UPGMA/WPGMA should recover it exactly
ULTRA = {"A": 2, "B": 2, "C": 4, "D": 4, "AB": 3, "CD": 1}
# Wikipedia "Neighbor joining" worked example: a2 b3 c4 d2 e1, internal 3 and 2
WIKI_NJ = {"a": 2, "b": 3, "c": 4, "d": 2, "e": 1, "ab": 3, "de": 2}
# Mount, Bioinformatics, Fitch-Margoliash worked example: A10 B12 C9 D9, internal 20
MOUNT_FM = {"A": 10, "B": 12, "C": 9, "D": 9, "AB": 20}


# --- the examples written in the scripts give the published answers ------------------------------

def test_upgma_matrix_example_is_wikipedia_wpgma(run_with):
    t = named(run_with("upgma_matrix.py")["tree"], "abcde")
    assert siblings(t, "a", "b") and siblings(t, "c", "d")
    assert leaf_lengths(t) == pytest.approx({"a": 8.5, "b": 8.5, "e": 11, "c": 14, "d": 14})
    assert all(t.distance(n) == pytest.approx(17.5) for n in "abcde")    # WPGMA root height (UPGMA: 16.5)


def test_nj_matrix_example_is_wikipedia(run_with):
    t = named(run_with("neighbor_joining_matrix.py")["tree"], "abcde")
    assert leaf_lengths(t) == pytest.approx({"a": 2, "b": 3, "c": 4, "d": 2, "e": 1})
    assert_paths(t, WIKI_NJ)


def test_fitch_matrix_example_is_mount(run_with):
    t = named(run_with("fitch_margoliash_matrix.py")["tree"], "ABCD")
    assert siblings(t, "A", "B") and siblings(t, "C", "D")
    assert leaf_lengths(t) == pytest.approx({"A": 10, "B": 12, "C": 9, "D": 9})
    assert_paths(t, MOUNT_FM)                                             # so internal branch = 20


@pytest.mark.parametrize("method", METHODS)
def test_seqs_example_runs(run_with, method):
    g = run_with(f"{method}_seqs.py")
    assert g["D"]["Seq1", "Seq2"] == 1                                    # 1 differing site
    assert siblings(g["tree"], "Seq1", "Seq2")


# --- matrix input --------------------------------------------------------------------------------

@pytest.mark.parametrize("method", METHODS)
def test_lower_triangle_same_as_full_matrix(run_with, method):
    full = [[0, 5, 9, 9, 8], [5, 0, 10, 10, 9], [9, 10, 0, 8, 7], [9, 10, 8, 0, 3], [8, 9, 7, 3, 0]]
    lower = [row[:i + 1] for i, row in enumerate(full)]
    a = run_with(f"{method}_matrix.py", matrix=full)["tree"]
    b = run_with(f"{method}_matrix.py", matrix=lower)["tree"]
    for x, y in combinations([f"Seq{i}" for i in range(1, 6)], 2):
        assert a.distance(x, y) == pytest.approx(b.distance(x, y))


@pytest.mark.parametrize("make", [from_seqs, from_matrix])
def test_upgma_ultrametric(run_with, make):
    t = make(run_with, "upgma", ULTRA)
    assert siblings(t, "A", "B") and siblings(t, "C", "D")
    assert leaf_lengths(t) == pytest.approx({"A": 2, "B": 2, "C": 4, "D": 4})
    assert all(t.distance(n) == pytest.approx(5) for n in "ABCD")        # all leaves at height 5
    assert_paths(t, ULTRA)


@pytest.mark.parametrize("make", [from_seqs, from_matrix])
def test_nj_recovers_additive_tree(run_with, make):
    t = make(run_with, "neighbor_joining", WIKI_NJ)
    assert leaf_lengths(t) == pytest.approx({n: k for n, k in WIKI_NJ.items() if len(n) == 1})
    assert_paths(t, WIKI_NJ)


@pytest.mark.parametrize("make", [from_seqs, from_matrix])
@pytest.mark.parametrize("edges", [MOUNT_FM, WIKI_NJ])
def test_fitch_recovers_tree(run_with, make, edges):
    assert_paths(make(run_with, "fitch_margoliash", edges), edges)


# --- sequence input ------------------------------------------------------------------------------

def test_distances_are_mismatch_counts(run_with):
    g = run_with("upgma_seqs.py", seqs=["ACGTACGTAC", "ACGTACGTTC", "TCGTACGTTA"])
    assert (g["D"]["Seq1", "Seq2"], g["D"]["Seq1", "Seq3"], g["D"]["Seq2", "Seq3"]) == (1, 3, 2)


def test_upgma_is_wpgma(run_with):
    # Biopython's upgma averages the 2 merged clusters equally (WPGMA, as taught in class).
    # merges: (A,B) h=.5 ; (AB,C) d=(3+2)/2=2.5 -> h=1.25 ; (ABC,D): WPGMA d=((10+9)/2+7)/2=8.25 -> h=4.125
    g = run_with("upgma_seqs.py", seqs=["AAAAAAAAAA", "AAAAAAAAAT", "AAAAAAATTT", "TTTTTTTTTT"])
    assert g["tree"].distance("Seq1") == pytest.approx(4.125)           # true UPGMA would give 4.333


@pytest.mark.parametrize("method", METHODS)
def test_auto_names(run_with, method):
    t = run_with(f"{method}_seqs.py", seqs=["ACGT", "ACGA", "TCGA"])["tree"]
    assert sorted(c.name for c in t.get_terminals()) == ["Seq1", "Seq2", "Seq3"]


# --- output --------------------------------------------------------------------------------------

@pytest.mark.parametrize("script", SCRIPTS)
def test_plot_shows_branch_lengths(run_with, script):
    plt.close("all")
    g = run_with(script)
    tree = g["tree"]
    labels = {t.get_text() for t in plt.gca().texts}
    for c in tree.find_clades():
        if round(c.branch_length or 0, 4):
            assert str(round(c.branch_length, 4)) in labels
    assert {c.name for c in tree.get_terminals()} <= {l.strip() for l in labels}   # leaves labelled
    if script.endswith("_seqs.py"):
        assert str(g["D"]) in g["_printed"]                               # distance matrix printed
