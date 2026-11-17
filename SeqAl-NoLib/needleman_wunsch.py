"""Needleman-Wunsch: GLOBAL alignment (all of both sequences is aligned).

No imports. Scores are whole numbers; `gap` is the cost of one gap position, so use a negative number.

    python needleman_wunsch.py
"""


def needleman_wunsch(a, b, match=1, mismatch=-1, gap=-2):
    """Returns (traceback, alignments).

    traceback[i][j] holds the arrows into cell (i, j), meaning "this cell's score came from ...":
        ↖  cell (i-1, j-1): a[i-1] is aligned with b[j-1]
        ↑  cell (i-1, j):   a[i-1] is aligned with a gap
        ←  cell (i, j-1):   b[j-1] is aligned with a gap
    A cell gets several arrows when routes tie. That is what produces several optimal alignments.

    alignments is a list of every optimal alignment as (row_a, row_b); "-" marks a gap.
    """
    n, m = len(a), len(b)

    # 1. Two tables with one extra row and column for "nothing aligned yet".
    #    Cell (i, j) is about the first i letters of a and the first j letters of b.
    scores = [[0] * (m + 1) for _ in range(n + 1)]
    traceback = [[""] * (m + 1) for _ in range(n + 1)]

    # 2. First column and first row: one sequence lined up against nothing but gaps.
    for i in range(1, n + 1):
        scores[i][0] = i * gap
        traceback[i][0] = "↑"
    for j in range(1, m + 1):
        scores[0][j] = j * gap
        traceback[0][j] = "←"

    # 3. Every other cell: the best of the three ways to arrive there.
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = scores[i-1][j-1] + (match if a[i-1] == b[j-1] else mismatch)
            up = scores[i-1][j] + gap
            left = scores[i][j-1] + gap
            best = max(diagonal, up, left)
            scores[i][j] = best
            if diagonal == best: traceback[i][j] += "↖"       # keep every route that ties for best
            if up == best:       traceback[i][j] += "↑"
            if left == best:     traceback[i][j] += "←"

    # 4. Walk back from the bottom-right corner to the top-left, following every arrow.
    #    We meet the alignment from its end, so each new letter goes on the FRONT of the rows.
    alignments = []

    def walk_back(i, j, row_a, row_b):
        if i == 0 and j == 0:                                  # reached the start: alignment is complete
            alignments.append((row_a, row_b))
        arrows = traceback[i][j]
        if "↖" in arrows: walk_back(i-1, j-1, a[i-1] + row_a, b[j-1] + row_b)
        if "↑" in arrows: walk_back(i-1, j,   a[i-1] + row_a, "-"    + row_b)
        if "←" in arrows: walk_back(i,   j-1, "-"    + row_a, b[j-1] + row_b)

    walk_back(n, m, "", "")
    return traceback, alignments


def show(a, b, traceback, alignments):
    print("Traceback matrix:   ↖ aligned pair   ↑ gap in b   ← gap in a   · start")
    print("    " + "".join(f"{letter:<4}" for letter in "-" + b))
    for letter, row in zip("-" + a, traceback):
        cells = [arrows or "·" for arrows in row]
        print(f"{letter}   " + "".join(f"{cell:<4}" for cell in cells))

    print(f"\n{len(alignments)} optimal alignment(s)")
    for number, (row_a, row_b) in enumerate(alignments, 1):
        bars = ""                                              # | marks each matching column
        for x, y in zip(row_a, row_b):
            bars += "|" if x == y else " "
        print(f"\n#{number}\n{row_a}\n{bars}\n{row_b}")


if __name__ == "__main__":
    a = input("Sequence 1: ").strip().upper()
    b = input("Sequence 2: ").strip().upper()
    match, mismatch, gap = map(int, input("Match, mismatch, gap scores (e.g. 1 -1 -2): ").split())
    traceback, alignments = needleman_wunsch(a, b, match, mismatch, gap)
    show(a, b, traceback, alignments)
