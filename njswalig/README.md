# njswalig

Needleman–Wunsch (global) and Smith–Waterman (local) pairwise sequence alignment in plain Python,
with no dependencies. It prints the alignment and, if you ask, the filled-in dynamic-programming score
table.

```bash
pip install njswalig
```

## Usage

```python
from njswalig import njalign      # Needleman-Wunsch (global)
from njswalig import swalign      # Smith-Waterman (local)

seqs = ["ATCG", "TTTA"]           # any number of sequences, any lengths
scores = [2, -1, -2]              # [match, mismatch, gap] -- always this order

alignment = njalign.sequences(seqs, scores)   # or njalign.sequence(...) -- same thing
alignment.show()
```

```
Needleman-Wunsch (global): Seq1 vs Seq2, score = -1
Seq1  ATCG
      .|..
Seq2  TTTA
```

`|` marks a match, `.` a mismatch and a space a gap.

### Smith–Waterman

```python
swalign.sequences(["TGTTACGG", "GGTTGACTA"], [3, -3, -2]).show()
```

```
Smith-Waterman (local): Seq1 vs Seq2, score = 13
Seq1  2 GTT-AC 6
        ||| ||
Seq2  2 GTTGAC 7
```

The numbers show where the aligned region starts and ends in each sequence (1-based).

### More than two sequences

Every pair is aligned. Sequences are named `Seq1`, `Seq2`, … in list order, or pass your own:

```python
njalign.sequences(["GATTACA", "GCATGCU", "GATCA"], [1, -1, -1], names=["a", "b", "c"]).show()
```

### The DP score table

```python
swalign.sequences(["ATCG", "TTTA"], [2, -1, -2]).show(matrix=True)
```

```
  - T T T A
- 0 0 0 0 0
A 0 0 0 0 2
T 0 2 2 2 0
C 0 0 1 1 1
G 0 0 0 0 0
```

### Using the results in code

`sequences(...)` returns a list of pair alignments. Each has `score`, `aligned1`, `aligned2`,
`match_line`, `name1`, `name2`, `seq1`, `seq2`, `matrix` (the DP table as a list of lists) and
1-based `start1`, `end1`, `start2`, `end2`.

```python
best = njalign.sequences(["GATTACA", "GCATGCU"], [1, -1, -1])[0]
print(best.score, best.aligned1, best.aligned2)     # 0 G-ATTACA GCA-TGCU
```

## Using your own matrix

Fill the score matrix yourself, then let njswalig check it, draw the traceback matrix and find
**every** optimal alignment:

```python
from njswalig import njalign, swalign

seq1, seq2 = "GATTACA", "GCATGCU"
match, mismatch, gap = 1, -1, -1

em = [[0] * (len(seq2) + 1) for i in range(len(seq1) + 1)]   # rows follow seq1, columns follow seq2
# ... fill em by hand or with your own loops ...

result = njalign.traceback(seq1, seq2, em, [match, mismatch, gap])   # swalign.traceback for local
result.show()
```

```
Needleman-Wunsch (global): Seq1 vs Seq2, match 1, mismatch -1, gap -1

Score matrix:
   -  G  C  A  T  G  C  U
-  0 -1 -2 -3 -4 -5 -6 -7
G -1  1  0 -1 -2 -3 -4 -5
...

Traceback matrix (↖ letter with letter, ↑ seq1 letter with a gap, ← seq2 letter with a gap):
   -  G  C  A  T  G  C  U
-  ·  ←  ←  ←  ←  ←  ←  ←
G  ↑  ↖  ←  ←  ← ↖←  ←  ←
...

3 optimal alignments, score 0:

#1  score = 0
Seq1  G-ATTACA
      | | |.|.
Seq2  GCA-TGCU
...
```

- `result.arrows`: the traceback matrix as a list of lists of strings (`""` where no arrow leaves a
  cell; ties give several arrows, in the order `↖↑←`).
- `result.alignments`: every optimal alignment, each with `score`, `aligned1`, `aligned2`,
  `match_line` and 1-based `start1`, `end1`, `start2`, `end2`.
- `result.score`: the optimal score.

**The matrix:** `em[i][j]` is the score for the first `i` letters of `seq1` against the first `j`
letters of `seq2`. So it has `len(seq1) + 1` rows of `len(seq2) + 1` numbers, border row and
column included. Global borders are `i * gap` and `j * gap`; local borders are 0, and no local cell
is below 0. A matrix built the other way round (transposed) is rejected with a message saying so.

**Checking:** every cell is recomputed from your scores. The first wrong one is named, for example
`cell (3, 2) is 4 but the scores give -1 (row 3 is seq1 letter 3 'T', column 2 is seq2 letter 2 'C')`.

**Many ties:** the number of optimal alignments can grow very fast (with a gap or mismatch score of
0, a short pair can have millions). `traceback(..., limit=100)` stops after 100; `result.complete`
then says whether that was all of them.

**Local convention:** the traceback starts from **every** cell holding the highest score and stops
at the first 0. This is the textbook method, and it keeps alignments that end in a run scoring 0 overall.
For example, with `AAA` vs `AACA` at `[1, -1, -1]`, both `AA`/`AA` and `AA-A`/`AACA` (match,
match, gap, match = 2) are reported. Biopython's `PairwiseAligner` leaves the second kind out, so
it can list fewer local alignments than njswalig.

## Details

- **Scoring** is linear: every gap position costs `gap`. `scores` must be a list of exactly three
  numbers, `[match, mismatch, gap]`; anything else raises `ValueError`.
- **Sequences** are stripped of surrounding spaces and upper-cased, so `"acgt "` equals `"ACGT"`.
- **Ties:** when several alignments share the best score, one is returned. The traceback prefers a
  diagonal step (match/mismatch), then a gap in the second sequence, then a gap in the first. Smith–Waterman
  starts from the first highest cell, scanning row by row.
- **No local match:** when nothing scores above 0, Smith–Waterman reports an empty alignment with score 0.

## License

MIT
