from Bio.Align import PairwiseAligner

def smith_waterman(a, b):     # LOCAL
    return PairwiseAligner(mode="local", match_score=1, mismatch_score=-1, gap_score=-2).align(a, b)[0]
