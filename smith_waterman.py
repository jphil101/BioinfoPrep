from Bio.Align import PairwiseAligner

a = input("Sequence 1: ").strip().upper()
b = input("Sequence 2: ").strip().upper()
match, mismatch, gap = map(float, input("Match, mismatch, gap scores (e.g. 1 -1 -2): ").split())
aligner = PairwiseAligner(mode="local", match_score=match, mismatch_score=mismatch, gap_score=gap)
best = aligner.align(a, b)[0]
print(best)
print("Score:", best.score)
