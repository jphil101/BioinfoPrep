"""Validation for both aligners. Plain asserts, no libraries.  Run: python test_seqal.py  (pytest works too)

`brute` enumerates every alignment by plain recursion (no DP, no traceback) as an independent oracle.
Local rule: an alignment starts with a match and no prefix scores <= 0 (else it would just start later).
"""
from needleman_wunsch import needleman_wunsch as nw
from smith_waterman import smith_waterman as sw

SCHEMES = [(1, -1, -2), (1, -1, -1), (2, -1, -1), (3, -2, -4), (1, 0, 0)]      # (match, mismatch, gap)


def brute(a, b, local, match, mismatch, gap):
    found = []                                                   # (score, start_i, start_j, row_a, row_b)

    def walk(si, sj, i, j, score, x, y):
        if score > 0 if local else (i, j) == (len(a), len(b)):
            found.append((score, si, sj, x, y))
        for di, dj in ((1, 1), (1, 0), (0, 1)):
            if i + di <= len(a) and j + dj <= len(b):
                s = score + (gap if not di * dj else match if a[i] == b[j] else mismatch)
                if not local or s > 0:
                    walk(si, sj, i + di, j + dj, s, x + (a[i] if di else "-"), y + (b[j] if dj else "-"))

    for si, sj in [(i, j) for i in range(len(a) + 1) for j in range(len(b) + 1)] if local else [(0, 0)]:
        walk(si, sj, si, sj, 0, "", "")
    top = max((f[0] for f in found), default=0)
    return sorted(f[1:] for f in found if f[0] == top)


def words(alphabet, max_len):                                    # every string of length 0..max_len
    out = [""]
    for w in out:                                                # list grows while we loop: breadth-first
        if len(w) < max_len: out += [w + c for c in alphabet]
    return out


def test_hand_checked():
    assert nw("A", "A") == ([["", "←"], ["↑", "↖"]], [("A", "A")])
    assert sorted(nw("AAA", "AA", 1, -1, -1)[1]) == [("AAA", "-AA"), ("AAA", "A-A"), ("AAA", "AA-")]   # 3-way tie
    assert sorted(nw("GATTACA", "GCATGCU", 1, -1, -1)[1]) == [                                  # textbook: score 0
        ("G-ATTACA", "GCA-TGCU"), ("G-ATTACA", "GCAT-GCU"), ("G-ATTACA", "GCATG-CU")]
    assert sorted(sw("ACAC", "AC")[1]) == [(0, 0, "AC", "AC"), (2, 0, "AC", "AC")]              # hit at two places


def test_edge_cases():
    assert nw("", "")[1] == [("", "")] and sw("", "")[1] == []
    assert nw("", "ACG")[1] == [("---", "ACG")] and nw("ACG", "")[1] == [("ACG", "---")]
    assert sw("", "ACG")[1] == [] and sw("ACG", "")[1] == []
    assert nw("ACGT", "ACGT")[1] == [("ACGT", "ACGT")] and sw("ACGT", "ACGT")[1] == [(0, 0, "ACGT", "ACGT")]
    assert sw("AAAA", "TTTT")[1] == []                           # nothing scores > 0
    assert len(sw("ACGT", "TGCA")[1]) == 4                       # four single-base hits tie
    t = nw("AC", "GT")[0]
    assert [len(r) for r in t] == [3, 3, 3] and t[0] == ["", "←", "←"] and [r[0] for r in t] == ["", "↑", "↑"]
    t = sw("AC", "GT")[0]
    assert t[0] == ["", "", ""] and [r[0] for r in t] == ["", "", ""]                           # local borders


def test_exhaustive_against_brute_force():
    for alphabet, length, schemes in (("AB", 4, SCHEMES), ("ACG", 3, SCHEMES[:1])):
        for scheme in schemes:
            for a in words(alphabet, length):
                for b in words(alphabet, length):
                    assert sorted((0, 0, x, y) for x, y in nw(a, b, *scheme)[1]) == brute(a, b, False, *scheme), (a, b, scheme)
                    assert sorted(sw(a, b, *scheme)[1]) == brute(a, b, True, *scheme), (a, b, scheme)


if __name__ == "__main__":
    for name, test in list(globals().items()):
        if name.startswith("test_"):
            test()
            print("ok  ", name)
