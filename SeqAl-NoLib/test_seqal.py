"""Checks for needleman_wunsch.py and smith_waterman.py.    Run:  python test_seqal.py    (pytest works too)

No libraries. The big test compares both aligners against `brute_force`, a deliberately naive
checker that tries EVERY possible alignment instead of filling a score table.
"""
from needleman_wunsch import needleman_wunsch
from smith_waterman import smith_waterman


def nw(a, b, *scores):                    # just the alignments, sorted so their order never matters
    return sorted(needleman_wunsch(a, b, *scores)[1])


def sw(a, b, *scores):
    return sorted(smith_waterman(a, b, *scores)[1])


def brute_force(a, b, local, match, mismatch, gap):
    """Try every way of lining up a and b; return the best ones as (start_i, start_j, row_a, row_b).
    Global: must use both sequences completely, from the start to the end.
    Local:  may start and stop anywhere, but never has a prefix scoring 0 or less (it would start later)."""
    found = []                                                  # (score, start_i, start_j, row_a, row_b)

    def extend(start_i, start_j, i, j, score, row_a, row_b):
        if local and row_a and score <= 0:
            return
        if (local and score > 0) or (not local and i == len(a) and j == len(b)):
            found.append((score, start_i, start_j, row_a, row_b))
        if i < len(a) and j < len(b):                           # pair a[i] with b[j]
            pair = match if a[i] == b[j] else mismatch
            extend(start_i, start_j, i + 1, j + 1, score + pair, row_a + a[i], row_b + b[j])
        if i < len(a):                                          # a[i] against a gap
            extend(start_i, start_j, i + 1, j, score + gap, row_a + a[i], row_b + "-")
        if j < len(b):                                          # gap against b[j]
            extend(start_i, start_j, i, j + 1, score + gap, row_a + "-", row_b + b[j])

    starts = [(0, 0)]
    if local:
        starts = [(i, j) for i in range(len(a) + 1) for j in range(len(b) + 1)]
    for i, j in starts:
        extend(i, j, i, j, 0, "", "")

    best = max((f[0] for f in found), default=0)
    return sorted(f[1:] for f in found if f[0] == best)


def all_strings(letters, max_length):                           # "", then every string up to max_length
    strings = [""]
    for s in strings:                                           # the list grows while we loop over it
        if len(s) < max_length:
            strings += [s + letter for letter in letters]
    return strings


# ---------------------------------------------------------------- known answers

def test_textbook_example():
    # The Wikipedia example: three alignments tie for the best score (0)
    assert nw("GATTACA", "GCATGCU", 1, -1, -1) == [
        ("G-ATTACA", "GCA-TGCU"),
        ("G-ATTACA", "GCAT-GCU"),
        ("G-ATTACA", "GCATG-CU")]


def test_every_tie_is_reported():
    # the spare A in "AAA" can be left unmatched in 3 different places
    assert nw("AAA", "AA", 1, -1, -1) == [("AAA", "-AA"), ("AAA", "A-A"), ("AAA", "AA-")]
    # "AC" occurs twice in "ACAC": same letters, different places
    assert sw("ACAC", "AC") == [(0, 0, "AC", "AC"), (2, 0, "AC", "AC")]
    # four different single-letter hits all score 1
    assert len(sw("ACGT", "TGCA")) == 4


def test_traceback_matrix():
    assert needleman_wunsch("A", "A")[0] == [["", "←"], ["↑", "↖"]]
    traceback = needleman_wunsch("AC", "GT")[0]
    assert len(traceback) == 3 and all(len(row) == 3 for row in traceback)      # (len(a)+1) x (len(b)+1)
    assert traceback[0] == ["", "←", "←"]                                       # global: first row is gaps
    assert [row[0] for row in traceback] == ["", "↑", "↑"]                      # global: first column is gaps
    traceback = smith_waterman("AC", "GT")[0]
    assert traceback == [["", "", ""], ["", "", ""], ["", "", ""]]              # local: nothing matches, no arrows


# ---------------------------------------------------------------- edge cases

def test_empty_sequences():
    assert nw("", "") == [("", "")]                             # global: the empty alignment
    assert nw("", "ACG") == [("---", "ACG")]                    # global: everything is a gap
    assert nw("ACG", "") == [("ACG", "---")]
    assert sw("", "") == sw("", "ACG") == sw("ACG", "") == []   # local: nothing to align


def test_identical_sequences():
    assert nw("ACGT", "ACGT") == [("ACGT", "ACGT")]
    assert sw("ACGT", "ACGT") == [(0, 0, "ACGT", "ACGT")]


def test_nothing_in_common():
    assert sw("AAAA", "TTTT") == []                             # local: no positive score, no alignment
    assert nw("AAAA", "TTTT") == [("AAAA", "TTTT")]             # global still has to align everything


def test_region_inside_longer_sequences():
    assert sw("TTACGTT", "GGACGGG") == [(2, 2, "ACG", "ACG")]   # only the shared ACG is reported


# ---------------------------------------------------------------- brute-force comparison

SCHEMES = [(1, -1, -2), (1, -1, -1), (2, -1, -1), (3, -2, -4), (1, 0, 0)]      # (match, mismatch, gap)


def compare_with_brute_force(strings, schemes):
    for scheme in schemes:
        for a in strings:
            for b in strings:
                global_alignments = [(0, 0, x, y) for x, y in nw(a, b, *scheme)]       # global always starts at (0, 0)
                assert global_alignments == brute_force(a, b, False, *scheme), ("NW", a, b, scheme)
                assert sw(a, b, *scheme) == brute_force(a, b, True, *scheme), ("SW", a, b, scheme)


def test_matches_brute_force_on_every_small_input():
    compare_with_brute_force(all_strings("AB", 4), SCHEMES)         # every pair of strings up to length 4
    compare_with_brute_force(all_strings("ACG", 3), SCHEMES[:1])    # a 3-letter alphabet, up to length 3


if __name__ == "__main__":
    for name, test in list(globals().items()):
        if name.startswith("test_"):
            test()
            print("ok  ", name)
