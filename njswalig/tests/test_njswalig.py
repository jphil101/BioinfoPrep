import random
from itertools import combinations

import pytest
from Bio.Align import PairwiseAligner

from njswalig import Aligner, Alignments, PairAlignment, njalign, swalign


def biopython_score(a, b, scores, mode):
    match, mismatch, gap = scores
    return PairwiseAligner(mode=mode, match_score=match, mismatch_score=mismatch, gap_score=gap).score(a, b)


def rescore(x, y, scores):
    match, mismatch, gap = scores
    return sum(gap if "-" in (p, q) else match if p == q else mismatch for p, q in zip(x, y))


random.seed(1)
CASES = [("".join(random.choices("ACGT", k=random.randint(1, 12))), "".join(random.choices("ACGT", k=random.randint(1, 12))),
          [random.randint(1, 4), random.randint(-4, 0), random.randint(-4, -1)]) for _ in range(150)]


# --- the usage from the brief ------------------------------------------------------------------------

def test_brief_usage(capsys):
    seqs = ["ATCG", "TTTA"]
    scores = [2, -1, -2]
    njalign.sequences(seqs, scores).show()
    option1 = capsys.readouterr().out
    njalign.sequence(seqs, scores).show()
    option2 = capsys.readouterr().out
    assert option1 == option2 == (
        "Needleman-Wunsch (global): Seq1 vs Seq2, score = -1\n"
        "Seq1  ATCG\n"
        "      .|..\n"
        "Seq2  TTTA\n")


def test_sequence_is_sequences():
    assert Aligner.sequence is Aligner.sequences


# --- correctness ---------------------------------------------------------------------------------------

@pytest.mark.parametrize("a,b,scores", CASES)
def test_needleman_wunsch_matches_biopython(a, b, scores):
    (pair,) = njalign.sequences([a, b], scores)
    assert pair.score == biopython_score(a, b, scores, "global")
    assert pair.aligned1.replace("-", "") == a and pair.aligned2.replace("-", "") == b   # whole sequences
    assert rescore(pair.aligned1, pair.aligned2, scores) == pair.score                    # alignment earns its score
    assert pair.matrix[-1][-1] == pair.score


@pytest.mark.parametrize("a,b,scores", CASES)
def test_smith_waterman_matches_biopython(a, b, scores):
    (pair,) = swalign.sequences([a, b], scores)
    assert pair.score == max(0, biopython_score(a, b, scores, "local"))
    assert pair.score == max(map(max, pair.matrix))
    if pair.score:
        assert a[pair.start1 - 1:pair.end1] == pair.aligned1.replace("-", "")            # region positions are right
        assert b[pair.start2 - 1:pair.end2] == pair.aligned2.replace("-", "")
        assert rescore(pair.aligned1, pair.aligned2, scores) == pair.score
    else:
        assert pair.aligned1 == pair.aligned2 == ""


def test_wikipedia_needleman_wunsch():
    (pair,) = njalign.sequences(["GATTACA", "GCATGCU"], [1, -1, -1])
    assert pair.score == 0


def test_wikipedia_smith_waterman():
    (pair,) = swalign.sequences(["TGTTACGG", "GGTTGACTA"], [3, -3, -2])
    assert (pair.aligned1, pair.aligned2, pair.score) == ("GTT-AC", "GTTGAC", 13)
    assert (pair.start1, pair.end1, pair.start2, pair.end2) == (2, 6, 2, 7)


def test_global_dp_table_edges():
    (pair,) = njalign.sequences(["AC", "GTA"], [1, -1, -2])
    assert [row[0] for row in pair.matrix] == [0, -2, -4]
    assert pair.matrix[0] == [0, -2, -4, -6]


def test_float_scores():
    (pair,) = njalign.sequences(["ACGT", "AGT"], [1.5, -0.5, -1.25])
    assert pair.score == pytest.approx(biopython_score("ACGT", "AGT", [1.5, -0.5, -1.25], "global"))


# --- many sequences ------------------------------------------------------------------------------------

@pytest.mark.parametrize("aligner", [njalign, swalign])
def test_every_pair_aligned(aligner):
    seqs = ["ACGT", "ACGA", "TCGA", "GGGT"]
    result = aligner.sequences(seqs, [2, -1, -2])
    assert isinstance(result, Alignments) and all(isinstance(p, PairAlignment) for p in result)
    assert [(p.name1, p.name2) for p in result] == [(f"Seq{i + 1}", f"Seq{j + 1}") for i, j in combinations(range(4), 2)]


def test_custom_names_and_cleanup():
    (pair,) = njalign.sequences(["  acgt ", "ACGA"], [1, -1, -1], names=["human", "mouse"])
    assert (pair.name1, pair.name2, pair.seq1) == ("human", "mouse", "ACGT")


def test_show_all_pairs_and_matrix(capsys):
    swalign.sequences(["ACG", "CG", "AG"], [2, -1, -2]).show(matrix=True)
    out = capsys.readouterr().out
    assert out.count("Smith-Waterman (local)") == 3
    assert "  - C G" in out                            # DP table header


def test_no_local_alignment_message(capsys):
    swalign.sequences(["AAAA", "TTTT"], [1, -1, -1]).show()
    assert "no positive-scoring local alignment" in capsys.readouterr().out


# --- input checking ------------------------------------------------------------------------------------

@pytest.mark.parametrize("scores", [[2, -1], [2, -1, -2, 0], "2 -1 -2", [2, "-1", -2], [True, -1, -2], None])
def test_scores_format_strictly_enforced(scores):
    with pytest.raises(ValueError, match=r"\[match, mismatch, gap\]"):
        njalign.sequences(["AC", "AG"], scores)


@pytest.mark.parametrize("seqs,msg", [("ACGT", "list of sequences"), (["ACGT"], "at least 2"),
                                      ([], "at least 2"), (["AC", 5], "list of sequences")])
def test_bad_sequences(seqs, msg):
    with pytest.raises(ValueError, match=msg):
        swalign.sequences(seqs, [1, -1, -1])


def test_bad_names_and_mode():
    with pytest.raises(ValueError, match="unique names"):
        njalign.sequences(["AC", "AG"], [1, -1, -1], names=["x", "x"])
    with pytest.raises(ValueError, match="global' or 'local"):
        Aligner("semi")
