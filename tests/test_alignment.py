import random
import pytest
from needleman_wunsch import needleman_wunsch
from smith_waterman import smith_waterman


def reference_score(a, b, local):
    """Plain DP (match 1, mismatch -1, gap -2) written independently of Biopython."""
    S = [[0 if local else -2 * (i + j) * (i == 0 or j == 0) for j in range(len(b) + 1)] for i in range(len(a) + 1)]
    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            v = max(S[i-1][j-1] + (1 if a[i-1] == b[j-1] else -1), S[i-1][j] - 2, S[i][j-1] - 2)
            S[i][j] = max(v, 0) if local else v
    return max(map(max, S)) if local else S[-1][-1]


def rescore(row1, row2):
    return sum(-2 if "-" in (x, y) else 1 if x == y else -1 for x, y in zip(row1, row2))


pairs = [("GATTACA", "GCATGCU"), ("ACGT", "ACGT"), ("AAAA", "TTTT"), ("ACACACTA", "AGCACACA")]
random.seed(0)
pairs += [("".join(random.choices("ACGT", k=random.randint(1, 12))),
           "".join(random.choices("ACGT", k=random.randint(1, 12)))) for _ in range(50)]


@pytest.mark.parametrize("a,b", pairs)
def test_needleman_wunsch(a, b):
    best = needleman_wunsch(a, b)
    assert best.score == reference_score(a, b, local=False)
    assert best[0].replace("-", "") == a and best[1].replace("-", "") == b   # global: whole sequences used
    assert rescore(best[0], best[1]) == best.score


@pytest.mark.parametrize("a,b", pairs)
def test_smith_waterman(a, b):
    if reference_score(a, b, local=True) == 0:      # no positive local match: nothing to align
        return
    best = smith_waterman(a, b)
    assert best.score == reference_score(a, b, local=True)
    assert best[0].replace("-", "") in a and best[1].replace("-", "") in b   # local: substrings
    assert rescore(best[0], best[1]) == best.score


def test_known_examples():
    assert needleman_wunsch("GATTACA", "GCATGCU").score == -1
    sw = smith_waterman("GATTACA", "GCATGCU")
    assert (sw[0], sw[1], sw.score) == ("AT", "AT", 2)
