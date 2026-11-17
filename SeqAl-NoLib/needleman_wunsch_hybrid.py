# Hybrid Needleman-Wunsch (global): the matrix is built by hand, then njswalig traces it back.

# ---------------------------------------------------------------- Part 1: the score matrix (plain Python)
seq1 = input("Enter the 1st sequence:\t").strip().upper()
seq2 = input("Enter the 2nd sequence:\t").strip().upper()
match = int(input("Enter the match score:\t"))
mismatch = int(input("Enter the mismatch score:\t"))
gap = int(input("Enter the gap penalty:\t"))

m = len(seq1)+1
n = len(seq2)+1


em = [[0]*n for i in range(m)]		# m rows (one per letter of seq1) x n columns (one per letter of seq2)

for i in range(m):
	em[i][0] = i*gap
for j in range(n):
	em[0][j] = j*gap

for i in range(1,m):
	for j in range(1,n):
		consider = []
		if seq1[i-1] == seq2[j-1]:
			consider.append(em[i-1][j-1]+match)
		if seq1[i-1] != seq2[j-1]:
			consider.append(em[i-1][j-1]+mismatch)
		consider.append(em[i-1][j]+gap)
		consider.append(em[i][j-1]+gap)
		em[i][j] = max(consider)


# ---------------------------------------------------------------- Part 2: njswalig takes over
# It checks every cell of em, then prints the score matrix, the traceback matrix and every optimal alignment.
try:
	from njswalig import njalign		# imported here, not at the top
	njalign.traceback
except (ImportError, AttributeError):
	raise SystemExit("\nNeeds njswalig 0.2.0 or later:  python3 -m pip install --upgrade njswalig")

result = njalign.traceback(seq1, seq2, em, [match, mismatch, gap])
result.show()
