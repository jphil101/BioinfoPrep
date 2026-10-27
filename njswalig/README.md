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
