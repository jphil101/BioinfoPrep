# Hybrid Smith-Waterman (local): the matrix is built by hand, then Biopython finds the alignments.

# ---------------------------------------------------------------- Part 1: the score matrix (plain Python)
seq1 = input("Enter the 1st sequence:\t")
seq2 = input("Enter the 2nd sequence:\t")
match = int(input("Enter the match score:\t"))
mismatch = int(input("Enter the mismatch score:\t"))
gap = int(input("Enter the gap penalty:\t"))

m = len(seq1)+1
n = len(seq2)+1


em = [[0]*n for i in range(m)]		# m rows (one per letter of seq1) x n columns (one per letter of seq2)
									# first row and column stay 0: no gap penalty before an alignment starts

for i in range(1,m):
	for j in range(1,n):
		consider = [0]				# the floor: a score never drops below 0
		if seq1[i-1] == seq2[j-1]:
			consider.append(em[i-1][j-1]+match)
		if seq1[i-1] != seq2[j-1]:
			consider.append(em[i-1][j-1]+mismatch)
		consider.append(em[i-1][j]+gap)
		consider.append(em[i][j-1]+gap)
		em[i][j] = max(consider)


# ---------------------------------------------------------------- Part 2: traceback arrows, read off the matrix
# An arrow is drawn when (neighbour's score + that move's cost) equals this cell's score; ties give several arrows.
# A cell scoring 0 gets no arrow: an alignment starts there.
arrows = [["·"]*n for i in range(m)]
for i in range(1,m):
	for j in range(1,n):
		cell = ""
		if em[i][j] > 0:
			if em[i][j] == em[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
				cell += "↖"
			if em[i][j] == em[i-1][j] + gap:
				cell += "↑"
			if em[i][j] == em[i][j-1] + gap:
				cell += "←"
		arrows[i][j] = cell or "·"


def show(title, matrix):
	width = max(len(str(cell)) for row in matrix for cell in row) + 1
	print(f"\n{title}")
	print("  " + "".join(f"{letter:>{width}}" for letter in "-" + seq2))
	for letter, row in zip("-" + seq1, matrix):
		print(f"{letter} " + "".join(f"{cell:>{width}}" for cell in row))


show("Score matrix:", em)
show("Traceback arrows:   ↖ aligned pair   ↑ gap in seq2   ← gap in seq1   · start", arrows)


# ---------------------------------------------------------------- Part 3: Biopython takes over
# Note: Biopython lists each local alignment in its shortest form. It leaves out a tied alignment that only
# extends a best alignment with a tail that nets 0 (e.g. a mismatch then a match), which the arrows above allow.
from Bio.Align import PairwiseAligner		# imported here, not at the top  (pip install biopython)

if not seq1 or not seq2:
	raise SystemExit("\nBiopython can't align an empty sequence.")

aligner = PairwiseAligner(mode="local", match_score=match, mismatch_score=mismatch, gap_score=gap)
alignments = aligner.align(seq1, seq2)
assert alignments.score == max(max(row) for row in em), "the hand-built matrix and Biopython disagree"	# cross-check

print(f"\n{len(alignments)} optimal alignment(s)")
for number, alignment in enumerate(alignments, 1):
	print(f"\n#{number}   score {int(alignment.score)}")
	print(alignment)
