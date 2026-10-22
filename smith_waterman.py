from Bio.Align import PairwiseAligner


# Smith-Waterman = LOCAL alignment
def smith_waterman(a, b, match=1, mismatch=-1, gap=-2):
    al = PairwiseAligner(mode="local", match_score=match, mismatch_score=mismatch, gap_score=gap)
    alns = al.align(a, b)
    return alns[0] if len(alns) else None          # None = nothing scores above 0


if __name__ == "__main__":
    best = smith_waterman("GATTACA", "GCATGCU")
    print(best, "score =", best.score)
