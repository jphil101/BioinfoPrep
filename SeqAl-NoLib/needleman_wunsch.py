"""Needleman-Wunsch: GLOBAL alignment (all of both sequences is aligned).

No imports. Scores are whole numbers; `gap` is the cost of one gap position, so use a negative number.

    python needleman_wunsch.py
"""


def needleman_wunsch(a, b, match=1, mismatch=-1, gap=-2):
    """Returns (scores, alignments).

    scores[i][j] is the best score for aligning the first i letters of a with the first j letters of b,
    so the bottom-right cell is the score of the whole alignment.
    alignments is a list of every optimal alignment as (row_a, row_b); "-" marks a gap.
    """
    n, m = len(a), len(b)

    def pair_score(i, j):                                      # lining up a[i-1] with b[j-1]
        return match if a[i-1] == b[j-1] else mismatch

    # 1. Table with one extra row and column for "nothing aligned yet".
    scores = [[0] * (m + 1) for _ in range(n + 1)]

    # 2. First column and first row: one sequence lined up against nothing but gaps.
    for i in range(1, n + 1): scores[i][0] = i * gap
    for j in range(1, m + 1): scores[0][j] = j * gap

    # 3. Every other cell: the best of the three ways to arrive there.
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = scores[i-1][j-1] + pair_score(i, j)     # a[i-1] with b[j-1]
            up = scores[i-1][j] + gap                          # a[i-1] with a gap
            left = scores[i][j-1] + gap                        # b[j-1] with a gap
            scores[i][j] = max(diagonal, up, left)

    # 4. Walk back from the bottom-right corner. A step back is valid if it explains this cell's score:
    #    the neighbour's score plus that move's cost equals it. If several steps are valid (a tie), follow
    #    them all: that is where multiple optimal alignments come from.
    #    We meet the alignment from its end, so each new letter goes on the FRONT of the rows.
    alignments = []

    def walk_back(i, j, row_a, row_b):
        if i == 0 and j == 0:                                  # reached the start: alignment is complete
            alignments.append((row_a, row_b))
        if i > 0 and j > 0 and scores[i][j] == scores[i-1][j-1] + pair_score(i, j):
            walk_back(i-1, j-1, a[i-1] + row_a, b[j-1] + row_b)
        if i > 0 and scores[i][j] == scores[i-1][j] + gap:
            walk_back(i-1, j,   a[i-1] + row_a, "-"    + row_b)
        if j > 0 and scores[i][j] == scores[i][j-1] + gap:
            walk_back(i,   j-1, "-"    + row_a, b[j-1] + row_b)

    walk_back(n, m, "", "")
    return scores, alignments


def show(a, b, scores, alignments):
    width = max(len(str(value)) for row in scores for value in row) + 1
    print("Score matrix:")
    print("  " + "".join(f"{letter:>{width}}" for letter in "-" + b))
    for letter, row in zip("-" + a, scores):
        print(f"{letter} " + "".join(f"{value:>{width}}" for value in row))

    score = scores[-1][-1]                                     # every optimal alignment has this score
    print(f"\n{len(alignments)} optimal alignment(s), score {score} each")
    for number, (row_a, row_b) in enumerate(alignments, 1):
        bars = ""                                              # | marks each matching column
        for x, y in zip(row_a, row_b):
            bars += "|" if x == y else " "
        print(f"\n#{number}   score {score}\n{row_a}\n{bars}\n{row_b}")


if __name__ == "__main__":
    a = input("Sequence 1: ").strip().upper()
    b = input("Sequence 2: ").strip().upper()
    match, mismatch, gap = map(int, input("Match, mismatch, gap scores (e.g. 1 -1 -2): ").split())
    scores, alignments = needleman_wunsch(a, b, match, mismatch, gap)
    show(a, b, scores, alignments)
