import pytest
import numpy as np


def blosum(run, *block):
    return run("blosum.py", str(len(block)), *block)["blosum"]


def pam(run, n, *seqs):
    return run("pam.py", str(len(seqs)), *seqs, str(n))["pam"]


def test_blosum_hand_calculated(run):
    # col1 AAA -> AA x3 ; col2 ABB -> AB x2, BB x1   => q: AA .5, AB 1/3, BB 1/6 ; p: A 2/3, B 1/3
    # s = 2 log2(q/e): AA 2log2(.5/(4/9)) = 0.34 ; AB 2log2((1/3)/(4/9)) = -0.83 ; BB 2log2((1/6)/(1/9)) = 1.17
    assert blosum(run, "AA", "AB", "AB").values.tolist() == [[0, -1], [-1, 1]]


def test_blosum_properties(run):
    S = blosum(run, "ACDA", "ACEA", "SCDA", "ACDG", "ACDA")
    assert (S.values == S.values.T).all()                     # symmetric
    assert list(S.index) == list(S.columns) == list("ACDEGS")
    assert S.loc["C", "C"] > 0                                # conserved column -> positive


def test_pam1_hand_calculated(run):
    # pairs (A,A),(A,B) both directions: q = [[1/2,1/4],[1/4,0]] ; f = [3/4, 1/4]
    # q/f = [[2/3,1/3],[1,0]] ; change = 1 - trace(q) = 1/2 ; lambda = 0.01/0.5 = 0.02
    # PAM1 = [[1-0.02/3, 0.02/3],[0.02, 0.98]] ; s = 10 log10(M_ij / f_j)
    #   AA 10log10(.99333/.75)=1.2  AB 10log10(.00667/.25)=-15.7  BA 10log10(.02/.75)=-15.7  BB 10log10(.98/.25)=5.9
    g = run("pam.py", "2", "AA", "AB", "1")
    assert np.allclose(g["M1"], [[1 - 0.02 / 3, 0.02 / 3], [0.02, 0.98]])
    assert g["pam"].values.tolist() == [[1, -16], [-16, 6]]


SEQS = ["ACDEAKLM", "ACDEAKLV", "ACEEAKLM", "SCDEAKIM", "ACDEGKLM", "ACDDAKLM"]


def test_pam1_properties(run):
    g = run("pam.py", str(len(SEQS)), *SEQS, "1")
    M1, f = np.asarray(g["M1"]), g["f"].values
    assert M1.sum(axis=1) == pytest.approx(1)            # rows are probabilities
    assert (f * (1 - np.diag(M1))).sum() == pytest.approx(0.01)   # exactly 1% change


def test_pam_scores(run):
    S = pam(run, 250, *SEQS)
    assert list(S.index) == list(S.columns) == sorted(set("".join(SEQS)))
    assert (S.values == S.values.T).all()                                # symmetric
    assert (np.diag(S.values) > 0).all()                                 # identities positive
    assert np.isfinite(S.loc["D", "E"]) and S.loc["D", "E"] > S.loc["D", "K"]   # observed swap beats never-seen


def test_pam_converges_to_zero(run):
    # PAM infinity: all info lost, M^n -> f, so every score -> 10 log10(1) = 0
    assert (pam(run, 5000, "AA", "AB").values == 0).all()
