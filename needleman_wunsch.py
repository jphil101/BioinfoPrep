from Bio.Align import PairwiseAligner


# Needleman-Wunsch = GLOBAL alignment
def needleman_wunsch(a, b, match=1, mismatch=-1, gap=-2):
    al = PairwiseAligner(mode="global", match_score=match, mismatch_score=mismatch, gap_score=gap)
    return al.align(a, b)[0]


if __name__ == "__main__":
    best = needleman_wunsch("GATTACA", "GCATGCU")
    print(best, "score =", best.score)
