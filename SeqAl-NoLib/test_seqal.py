"""Validation for seqal.py. Plain asserts, only import is the module under test.

    python test_seqal.py        (pytest picks up the test_ functions too)

The oracle `brute` enumerates every alignment by plain recursion (no DP, no traceback), so it
checks the DP + traceback independently. Local rule matches Smith-Waterman: every alignment starts
with a match and no prefix scores <= 0 (otherwise it would simply start later).
"""
from seqal import align

SCHEMES = [(1, -1, -2), (1, -1, -1), (2, -1, -1), (3, -2, -4), (1, 0, 0)]   # (match, mismatch, gap)


def brute(a, b, local, match, mismatch, gap):
    found = []                                            # (score, start_i, start_j, row_a, row_b)

    def walk(si, sj, i, j, score, x, y):
        if score > 0 if local else (i, j) == (len(a), len(b)):
            found.append((score, si, sj, x, y))
        for di, dj in ((1, 1), (1, 0), (0, 1)):
            if i + di > len(a) or j + dj > len(b): continue
            s = score + gap if not (di and dj) else score + (match if a[i] == b[j] else mismatch)
            if local and s <= 0: continue
            walk(si, sj, i + di, j + dj, s, x + (a[i] if di else "-"), y + (b[j] if dj else "-"))

    for si, sj in [(i, j) for i in range(len(a) + 1) for j in range(len(b) + 1)] if local else [(0, 0)]:
        walk(si, sj, si, sj, 0, "", "")
    top = max((f[0] for f in found), default=0)
    return sorted(f[1:] for f in found if f[0] == top)


def words(alphabet, max_len):
    out = [""]
    for w in out:                                         # list grows while we loop: breadth-first
        if len(w) < max_len: out += [w + c for c in alphabet]
    return out


def test_hand_checked():
    assert align("A", "A") == ([["", "L"], ["U", "D"]], [(0, 0, "A", "A")])
    # 3 gap positions tie (2 matches, 1 gap)
    assert sorted(align("AAA", "AA", False, 1, -1, -1)[1]) == [(0, 0, "AAA", "-AA"), (0, 0, "AAA", "A-A"), (0, 0, "AAA", "AA-")]
    # same local hit at two places in a
    assert sorted(align("ACAC", "AC", True)[1]) == [(0, 0, "AC", "AC"), (2, 0, "AC", "AC")]
    # textbook example (Wikipedia): score 0, three optimal alignments
    got = align("GATTACA", "GCATGCU", False, 1, -1, -1)[1]
    assert sorted(g[2:] for g in got) == [("G-ATTACA", "GCA-TGCU"), ("G-ATTACA", "GCAT-GCU"), ("G-ATTACA", "GCATG-CU")]


def test_edge_cases():
    assert align("", "")[1] == [(0, 0, "", "")]           # global: empty alignment
    assert align("", "", True)[1] == []                   # local: nothing to align
    assert align("", "ACG")[1] == [(0, 0, "---", "ACG")]  # one side empty: all gaps
    assert align("ACG", "")[1] == [(0, 0, "ACG", "---")]
    assert align("", "ACG", True)[1] == []
    assert align("ACGT", "ACGT")[1] == [(0, 0, "ACGT", "ACGT")]
    assert align("AAAA", "TTTT", True)[1] == []           # nothing scores > 0
    assert len(align("ACGT", "TGCA", True)[1]) == 4       # four separate single-base hits tie
    trace = align("AC", "GT")[0]
    assert [len(r) for r in trace] == [3, 3, 3]           # (len(a)+1) x (len(b)+1)
    assert trace[0] == ["", "L", "L"] and [r[0] for r in trace] == ["", "U", "U"]   # global borders
    trace = align("AC", "GT", True)[0]
    assert trace[0] == ["", "", ""] and [r[0] for r in trace] == ["", "", ""]       # local borders: alignment starts


def test_exhaustive_against_brute_force():
    seqs = words("AB", 4)                                 # every string of length 0..4 over {A, B}
    for scheme in SCHEMES:
        for local in (False, True):
            for a in seqs:
                for b in seqs:
                    got = sorted(align(a, b, local, *scheme)[1])
                    assert got == brute(a, b, local, *scheme), (a, b, local, scheme)
    seqs = words("ACG", 3)                                # three-letter alphabet, shorter
    for local in (False, True):
        for a in seqs:
            for b in seqs:
                assert sorted(align(a, b, local)[1]) == brute(a, b, local, 1, -1, -2), (a, b, local)


if __name__ == "__main__":
    for name, test in list(globals().items()):
        if name.startswith("test_"):
            test()
            print("ok  ", name)
