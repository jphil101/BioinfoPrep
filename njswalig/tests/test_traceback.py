"""Tests for Aligner.traceback: hand-built score matrix -> traceback matrix + every optimal alignment."""
import pytest

from njswalig import PairAlignment, Traceback, njalign, swalign


# ---------------------------------------------------------------- helpers (independent of the library)

def build_matrix(a, b, local, match, mismatch, gap):
    """Fill the score matrix the way a student does by hand (rows follow a, columns follow b)."""
    rows, cols = len(a) + 1, len(b) + 1
    em = [[0] * cols for _ in range(rows)]
    if not local:
        for i in range(rows):
            em[i][0] = i * gap
        for j in range(cols):
            em[0][j] = j * gap
    for i in range(1, rows):
        for j in range(1, cols):
            pair = match if a[i - 1] == b[j - 1] else mismatch
            em[i][j] = max(em[i - 1][j - 1] + pair, em[i - 1][j] + gap, em[i][j - 1] + gap)
            if local:
                em[i][j] = max(em[i][j], 0)
    return em


def trace(a, b, local, match, mismatch, gap, **kw):
    aligner = swalign if local else njalign
    return aligner.traceback(a, b, build_matrix(a, b, local, match, mismatch, gap), [match, mismatch, gap], **kw)


def brute_force(a, b, local, match, mismatch, gap):
    """Try every way of lining up a and b -- no score table, no traceback.
    Returns (best score, sorted best alignments as (start_i, start_j, row_a, row_b)).
    Global: uses both sequences completely. Local: may start and stop anywhere, but never has a
    prefix scoring 0 or less (it would start later instead)."""
    found = []

    def extend(start_i, start_j, i, j, score, row_a, row_b):
        if local and row_a and score <= 0:
            return
        if (local and score > 0) or (not local and i == len(a) and j == len(b)):
            found.append((score, start_i, start_j, row_a, row_b))
        if i < len(a) and j < len(b):
            pair = match if a[i] == b[j] else mismatch
            extend(start_i, start_j, i + 1, j + 1, score + pair, row_a + a[i], row_b + b[j])
        if i < len(a):
            extend(start_i, start_j, i + 1, j, score + gap, row_a + a[i], row_b + "-")
        if j < len(b):
            extend(start_i, start_j, i, j + 1, score + gap, row_a + "-", row_b + b[j])

    starts = [(i, j) for i in range(len(a) + 1) for j in range(len(b) + 1)] if local else [(0, 0)]
    for i, j in starts:
        extend(i, j, i, j, 0, "", "")
    best = max((f[0] for f in found), default=0)
    return best, sorted(f[1:] for f in found if f[0] == best)


def all_strings(letters, max_length):
    strings = [""]
    for s in strings:
        if len(s) < max_length:
            strings += [s + letter for letter in letters]
    return strings


def as_tuples(result):
    """(start_i, start_j, aligned1, aligned2), 0-based starts, sorted -- duplicates stay visible."""
    return sorted((p.start1 - 1, p.start2 - 1, p.aligned1, p.aligned2) for p in result.alignments)


def pairs(result):
    return sorted((p.aligned1, p.aligned2) for p in result.alignments)


SCHEMES = [(1, -1, -2), (1, -1, -1), (2, -1, -1), (3, -2, -4), (1, 0, 0)]      # (match, mismatch, gap)


# ---------------------------------------------------------------- 1. exhaustive vs the brute-force oracle

@pytest.mark.parametrize("scheme", SCHEMES)
@pytest.mark.parametrize("local", [False, True], ids=["global", "local"])
def test_matches_brute_force_on_every_small_input(scheme, local):
    strings = all_strings("AB", 4)                              # "", "A", "B", "AA", ... "BBBB"
    for a in strings:
        for b in strings:
            result = trace(a, b, local, *scheme)
            best, expected = brute_force(a, b, local, *scheme)
            where = (a, b, scheme, "local" if local else "global")
            assert as_tuples(result) == expected, where
            assert all(p.score == best for p in result.alignments), where
            if expected:
                assert result.score == best, where


def test_matches_brute_force_three_letters():
    for a in all_strings("ACG", 3):
        for b in all_strings("ACG", 3):
            for local in (False, True):
                assert as_tuples(trace(a, b, local, 1, -1, -2)) == brute_force(a, b, local, 1, -1, -2)[1]


# ---------------------------------------------------------------- 2. edge cases

def test_both_empty():
    result = trace("", "", False, 1, -1, -2)
    assert pairs(result) == [("", "")] and result.alignments[0].score == 0
    assert trace("", "", True, 1, -1, -2).alignments == []


def test_one_side_empty():
    assert pairs(trace("", "ACG", False, 1, -1, -2)) == [("---", "ACG")]
    assert pairs(trace("ACG", "", False, 1, -1, -2)) == [("ACG", "---")]
    assert trace("ACG", "", True, 1, -1, -2).alignments == []
    assert trace("", "ACG", True, 1, -1, -2).alignments == []


def test_identical_sequences():
    assert pairs(trace("ACGT", "ACGT", False, 1, -1, -2)) == [("ACGT", "ACGT")]
    assert as_tuples(trace("ACGT", "ACGT", True, 1, -1, -2)) == [(0, 0, "ACGT", "ACGT")]


def test_local_nothing_in_common():
    result = trace("AAAA", "TTTT", True, 1, -1, -1)
    assert result.alignments == [] and result.score == 0
    assert all(cell == "" for row in result.arrows for cell in row)


def test_every_tie_is_reported():
    result = trace("AAA", "AA", False, 1, -1, -1)
    assert pairs(result) == [("AAA", "-AA"), ("AAA", "A-A"), ("AAA", "AA-")]
    assert [p.score for p in result.alignments] == [1, 1, 1]


def test_local_same_letters_different_places():
    assert as_tuples(trace("ACAC", "AC", True, 1, -1, -2)) == [(0, 0, "AC", "AC"), (2, 0, "AC", "AC")]


def test_textbook_gattaca():
    result = trace("GATTACA", "GCATGCU", False, 1, -1, -1)
    assert result.score == 0
    assert pairs(result) == [("G-ATTACA", "GCA-TGCU"), ("G-ATTACA", "GCAT-GCU"), ("G-ATTACA", "GCATG-CU")]


def test_unequal_lengths():
    result = trace("ACGTAC", "AGC", False, 2, -1, -2)
    assert len(result.arrows) == 7 and all(len(row) == 4 for row in result.arrows)
    for p in result.alignments:
        assert p.aligned1.replace("-", "") == "ACGTAC" and p.aligned2.replace("-", "") == "AGC"


def test_zero_net_tail_kept():
    # match, match, gap, match scores the same as match, match: the textbook traceback keeps both
    assert ("AA-A", "AACA") in pairs(trace("AAA", "AACA", True, 1, -1, -1))


def test_local_positions():
    result = trace("TTACGTT", "GGACGGG", True, 1, -1, -2)
    (p,) = result.alignments
    assert (p.aligned1, p.start1, p.end1, p.start2, p.end2) == ("ACG", 3, 5, 3, 5)


def test_arrows_matrix():
    result = trace("AC", "AG", False, 1, -1, -2)
    assert result.arrows == [["", "←", "←"],
                             ["↑", "↖", "←"],
                             ["↑", "↑", "↖"]]
    assert result.matrix == [[0, -2, -4], [-2, 1, -1], [-4, -1, 0]]


def test_case_and_spaces_normalised():
    em = build_matrix("ACG", "AG", False, 1, -1, -2)
    assert pairs(njalign.traceback("  acg ", "ag", em, [1, -1, -2])) == pairs(njalign.traceback("ACG", "AG", em, [1, -1, -2]))


def test_float_scores():
    result = trace("ACGT", "AGT", False, 1.5, -0.5, -1.25)
    assert result.alignments and all(p.score == pytest.approx(result.score) for p in result.alignments)


def test_result_types_and_names():
    result = njalign.traceback("AC", "AG", build_matrix("AC", "AG", False, 1, -1, -2), [1, -1, -2], names=["x", "y"])
    assert isinstance(result, Traceback) and all(isinstance(p, PairAlignment) for p in result.alignments)
    assert (result.name1, result.name2, result.alignments[0].name1) == ("x", "y", "x")


def test_limit_and_exponential_ties():
    # nothing matches and gap 0 = mismatch 0: every path through the table ties (8,097,453 of them)
    full = trace("A" * 10, "B" * 10, False, 1, 0, 0, limit=50)
    assert len(full.alignments) == 50 and not full.complete
    assert "there are more" in str(full)
    assert trace("AAA", "AA", False, 1, -1, -1, limit=3).complete
    assert trace("AAA", "AA", False, 1, -1, -1, limit=10).complete
    assert len(trace("AAA", "AA", False, 1, -1, -1, limit=2).alignments) == 2


def test_long_sequences_no_recursion_error():
    a = "ACGT" * 500                                            # 2000 letters
    b = "ACGA" * 500
    for local in (False, True):
        result = trace(a, b, local, 1, -1, -2, limit=1)
        assert result.alignments and result.alignments[0].score == result.score


def test_sequences_unchanged():
    # the existing API still returns one alignment per pair
    (pair,) = njalign.sequences(["GATTACA", "GCATGCU"], [1, -1, -1])
    assert pair.score == 0


# ---------------------------------------------------------------- 3. contract errors

def test_wrong_shape():
    with pytest.raises(ValueError, match="must have 4 rows"):
        njalign.traceback("ACG", "AG", [[0, -2, -4], [-2, 1, -1]], [1, -1, -2])
    with pytest.raises(ValueError, match="must have 4 rows"):
        njalign.traceback("ACG", "AG", [[0, -2], [-2, 1], [-4, -1], [-6, -3]], [1, -1, -2])


def test_transposed_matrix():
    em = build_matrix("ACG", "AG", False, 1, -1, -2)
    transposed = [list(column) for column in zip(*em)]          # 3 rows of 4: the classic mistake
    with pytest.raises(ValueError, match="transposed"):
        njalign.traceback("ACG", "AG", transposed, [1, -1, -2])


@pytest.mark.parametrize("local", [False, True])
def test_tampered_cell_is_named(local):
    em = build_matrix("GATTACA", "GCATGCU", local, 1, -1, -1)
    em[3][2] += 3
    aligner = swalign if local else njalign
    with pytest.raises(ValueError, match=r"cell \(3, 2\) is -?\d+ but the scores give -?\d+"):
        aligner.traceback("GATTACA", "GCATGCU", em, [1, -1, -1])


def test_wrong_border():
    em = build_matrix("AC", "AG", False, 1, -1, -2)
    em[2][0] = 4                                                # border should be 2 * gap = -4
    with pytest.raises(ValueError, match=r"cell \(2, 0\) is 4 but the scores give -4"):
        njalign.traceback("AC", "AG", em, [1, -1, -2])
    em = build_matrix("AC", "AG", True, 1, -1, -2)
    em[0][1] = -2                                               # local border must be 0
    with pytest.raises(ValueError, match=r"cell \(0, 1\) is -2 but the scores give 0"):
        swalign.traceback("AC", "AG", em, [1, -1, -2])


def test_global_matrix_given_to_local():
    em = build_matrix("AC", "AG", False, 1, -1, -2)             # has negative borders
    with pytest.raises(ValueError, match=r"cell \(0, 1\)"):
        swalign.traceback("AC", "AG", em, [1, -1, -2])


@pytest.mark.parametrize("scores", [[1, -1], [1, -1, -2, 0], "1 -1 -2", [1, "-1", -2], [True, -1, -2]])
def test_bad_scores(scores):
    with pytest.raises(ValueError, match=r"\[match, mismatch, gap\]"):
        njalign.traceback("AC", "AG", build_matrix("AC", "AG", False, 1, -1, -2), scores)


@pytest.mark.parametrize("matrix,msg", [("nope", "list of rows"), ([[0, "x"]], "only numbers"),
                                        ([[0, True]], "only numbers")])
def test_bad_matrix_type(matrix, msg):
    with pytest.raises(ValueError, match=msg):
        njalign.traceback("A", "", matrix, [1, -1, -2])


@pytest.mark.parametrize("limit", [0, -1, 1.5, True, "3"])
def test_bad_limit(limit):
    with pytest.raises(ValueError, match="limit"):
        njalign.traceback("AC", "AG", build_matrix("AC", "AG", False, 1, -1, -2), [1, -1, -2], limit=limit)


def test_bad_sequences():
    with pytest.raises(ValueError, match="must be strings"):
        njalign.traceback(["AC"], "AG", [[0]], [1, -1, -2])


# ---------------------------------------------------------------- 4. exact printed output

def test_show_global(capsys):
    trace("GATTACA", "GCATGCU", False, 1, -1, -1).show()
    assert capsys.readouterr().out == """\
Needleman-Wunsch (global): Seq1 vs Seq2, match 1, mismatch -1, gap -1

Score matrix:
   -  G  C  A  T  G  C  U
-  0 -1 -2 -3 -4 -5 -6 -7
G -1  1  0 -1 -2 -3 -4 -5
A -2  0  0  1  0 -1 -2 -3
T -3 -1 -1  0  2  1  0 -1
T -4 -2 -2 -1  1  1  0 -1
A -5 -3 -3 -1  0  0  0 -1
C -6 -4 -2 -2 -1 -1  1  0
A -7 -5 -3 -1 -2 -2  0  0

Traceback matrix (↖ letter with letter, ↑ seq1 letter with a gap, ← seq2 letter with a gap):
   -  G  C  A  T  G  C  U
-  ·  ←  ←  ←  ←  ←  ←  ←
G  ↑  ↖  ←  ←  ← ↖←  ←  ←
A  ↑  ↑  ↖  ↖  ←  ←  ←  ←
T  ↑  ↑ ↖↑  ↑  ↖  ←  ←  ←
T  ↑  ↑ ↖↑  ↑ ↖↑  ↖ ↖← ↖←
A  ↑  ↑ ↖↑  ↖  ↑ ↖↑  ↖ ↖←
C  ↑  ↑  ↖  ↑  ↑ ↖↑  ↖  ←
A  ↑  ↑  ↑  ↖ ↑← ↖↑  ↑  ↖

3 optimal alignments, score 0:

#1  score = 0
Seq1  G-ATTACA
      | | |.|.
Seq2  GCA-TGCU

#2  score = 0
Seq1  G-ATTACA
      | || .|.
Seq2  GCAT-GCU

#3  score = 0
Seq1  G-ATTACA
      | ||. |.
Seq2  GCATG-CU
"""


def test_show_local(capsys):
    trace("ATCG", "TTTA", True, 2, -1, -2).show()
    assert capsys.readouterr().out == """\
Smith-Waterman (local): Seq1 vs Seq2, match 2, mismatch -1, gap -2

Score matrix:
  - T T T A
- 0 0 0 0 0
A 0 0 0 0 2
T 0 2 2 2 0
C 0 0 1 1 1
G 0 0 0 0 0

Traceback matrix (↖ letter with letter, ↑ seq1 letter with a gap, ← seq2 letter with a gap):
  - T T T A
- · · · · ·
A · · · · ↖
T · ↖ ↖ ↖ ·
C · · ↖ ↖ ↖
G · · · · ·

4 optimal alignments, score 2:

#1  score = 2
Seq1  1 A 1
        |
Seq2  4 A 4

#2  score = 2
Seq1  2 T 2
        |
Seq2  1 T 1

#3  score = 2
Seq1  2 T 2
        |
Seq2  2 T 2

#4  score = 2
Seq1  2 T 2
        |
Seq2  3 T 3
"""


def test_show_no_local_alignment(capsys):
    trace("AA", "TT", True, 1, -1, -1).show()
    assert capsys.readouterr().out.endswith("No local alignment: no cell scores above 0.\n")


def test_walk_returns_at_zero_even_if_an_arrow_points_on():
    # The walk itself must stop at a 0 cell (not only because 0 cells get no arrows): give it a
    # traceback matrix where a 0 cell does have an arrow, and the result must not change.
    from njswalig import _arrows, _walk
    a, b = "AA", "A"
    S = build_matrix(a, b, True, 1, -1, 0)                       # gap 0: [[0, 0], [0, 1], [0, 1]]
    arrows = _arrows("local", a, b, S, 1, -1, 0)
    assert arrows[1][0] == ""
    expected = _walk("local", a, b, S, arrows, None)
    arrows[1][0] = "↑"                                          # a 0 cell that points further up
    assert _walk("local", a, b, S, arrows, None) == expected
