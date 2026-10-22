import numpy as np
import pytest
from blosum import blosum
from pam import pam1, pam


def test_blosum_hand_calculated():
    # col1 AAA -> AA x3 ; col2 ABB -> AB x2, BB x1   => q: AA .5, AB 1/3, BB 1/6 ; p: A 2/3, B 1/3
    # s = 2 log2(q/e): AA 2log2(.5/(4/9)) = 0.34 ; AB 2log2((1/3)/(4/9)) = -0.83 ; BB 2log2((1/6)/(1/9)) = 1.17
    S = blosum(["AA", "AB", "AB"])
    assert S.values.tolist() == [[0, -1], [-1, 1]]


def test_blosum_properties():
    S = blosum(["ACDA", "ACEA", "SCDA", "ACDG", "ACDA"])
    assert (S.values == S.values.T).all()                     # symmetric
    assert list(S.index) == list(S.columns) == list("ACDEGS")
    assert S.loc["C", "C"] > 0                                # conserved column -> positive


def test_pam1_hand_calculated():
    # pairs (A,A),(A,B) both directions: F = [[2,1],[1,0]] ; f = [3/4, 1/4]
    # M = [[2/3,1/3],[1,0]] ; change = 3/4*1/3 + 1/4*1 = 1/2 ; lambda = 0.01/0.5 = 0.02
    M1, f = pam1(["AA", "AB"])
    assert f.tolist() == pytest.approx([0.75, 0.25])
    assert np.allclose(M1.values, [[1 - 0.02 / 3, 0.02 / 3], [0.02, 0.98]])


SEQS = ["ACDEAKLM", "ACDEAKLV", "ACEEAKLM", "SCDEAKIM", "ACDEGKLM", "ACDDAKLM"]


def test_pam1_properties():
    M1, f = pam1(SEQS)
    assert M1.sum(axis=1).values == pytest.approx(1)                    # rows are probabilities
    assert (f * (1 - np.diag(M1))).sum() == pytest.approx(0.01)         # exactly 1% change = 1 PAM
    assert f.values @ M1.values == pytest.approx(f.values)              # frequencies are stationary


def test_pam_scores():
    S = pam(SEQS, 250)
    finite = S.values[np.isfinite(S.values)]
    assert (S.values == S.values.T).all()                                # symmetric
    assert (np.diag(S.values) > 0).all()                                 # identities positive
    assert np.isfinite(S.loc["D", "E"]) and S.loc["D", "E"] > S.loc["D", "K"]   # observed swap beats never-seen
    assert len(finite) > 0


def test_pam_converges_to_zero():
    # PAM infinity: all info lost, M^n -> f, so every score -> 10 log10(1) = 0
    assert (pam(["AA", "AB"], 5000).values == 0).all()
