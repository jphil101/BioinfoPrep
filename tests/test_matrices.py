import numpy as np
import pytest
from blosum import blosum
from pam import pam


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
    # pairs (A,A),(A,B) both directions: q = [[1/2,1/4],[1/4,0]] ; f = [3/4, 1/4]
    # M = q/f = [[2/3,1/3],[1,0]] ; change = 1 - trace(q) = 1/2 ; lambda = 0.01/0.5 = 0.02
    # PAM1 = [[1-0.02/3, 0.02/3],[0.02, 0.98]] ; s = 10 log10(M_ij / f_j)
    #   AA 10log10(.99333/.75)=1.2  AB 10log10(.00667/.25)=-15.7  BA 10log10(.02/.75)=-15.7  BB 10log10(.98/.25)=5.9
    assert pam(["AA", "AB"], 1).values.tolist() == [[1, -16], [-16, 6]]


SEQS = ["ACDEAKLM", "ACDEAKLV", "ACEEAKLM", "SCDEAKIM", "ACDEGKLM", "ACDDAKLM"]


def test_pam_scores():
    S = pam(SEQS, 250)
    assert list(S.index) == list(S.columns) == sorted(set("".join(SEQS)))
    assert (S.values == S.values.T).all()                                # symmetric
    assert (np.diag(S.values) > 0).all()                                 # identities positive
    assert np.isfinite(S.loc["D", "E"]) and S.loc["D", "E"] > S.loc["D", "K"]   # observed swap beats never-seen


def test_pam0_is_identity():
    # PAM0: M = I, so diagonal = 10 log10(1/f_i), off-diagonal = -inf
    S = pam(["AA", "AB"], 0)
    assert S.loc["A", "A"] == round(10 * np.log10(1 / 0.75)) and S.loc["B", "B"] == round(10 * np.log10(4))


def test_pam_converges_to_zero():
    # PAM infinity: all info lost, M^n -> f, so every score -> 10 log10(1) = 0
    assert (pam(["AA", "AB"], 5000).values == 0).all()
