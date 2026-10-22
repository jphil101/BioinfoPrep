from Bio.Align import PairwiseAligner

def needleman_wunsch(a, b):   # GLOBAL
    return PairwiseAligner(mode="global", match_score=1, mismatch_score=-1, gap_score=-2).align(a, b)[0]
