"""Smith-Waterman: LOCAL alignment (finds the best-matching region; the rest of each sequence is ignored).

No imports. Scores are whole numbers; `gap` is the cost of one gap position, so use a negative number.

    python smith_waterman.py
"""


def smith_waterman(a, b, match=1, mismatch=-1, gap=-2):
    """Returns (traceback, alignments).

    traceback[i][j] holds the arrows into cell (i, j), meaning "this cell's score came from ...":
        ↖  cell (i-1, j-1): a[i-1] is aligned with b[j-1]
        ↑  cell (i-1, j):   a[i-1] is aligned with a gap
        ←  cell (i, j-1):   b[j-1] is aligned with a gap
    A cell gets several arrows when routes tie. That is what produces several optimal alignments.
    An EMPTY cell (score 0) has no arrow: an alignment starts there.

    alignments is a list of every optimal alignment as (start_i, start_j, row_a, row_b).
    It begins at a[start_i] and b[start_j]; "-" marks a gap.
    """
    n, m = len(a), len(b)

    # 1. Two tables with one extra row and column for "nothing aligned yet".
    #    Cell (i, j) is about the first i letters of a and the first j letters of b.
    #    Unlike Needleman-Wunsch, the first row and column stay 0: skipping the start of a sequence is free.
    scores = [[0] * (m + 1) for _ in range(n + 1)]
    traceback = [[""] * (m + 1) for _ in range(n + 1)]

    # 2. Every other cell: the best of the three ways to arrive there, but never below 0.
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = scores[i-1][j-1] + (match if a[i-1] == b[j-1] else mismatch)
            up = scores[i-1][j] + gap
            left = scores[i][j-1] + gap
            best = max(diagonal, up, left, 0)                  # 0 means "start a fresh alignment here"
            scores[i][j] = best
            if best > 0:                                       # (at 0 there is no arrow: nothing to trace back)
                if diagonal == best: traceback[i][j] += "↖"    # keep every route that ties for best
                if up == best:       traceback[i][j] += "↑"
                if left == best:     traceback[i][j] += "←"

    # 3. Walk back from every cell holding the highest score, until we hit a cell with no arrow.
    #    We meet the alignment from its end, so each new letter goes on the FRONT of the rows.
    alignments = []

    def walk_back(i, j, row_a, row_b):
        arrows = traceback[i][j]
        if arrows == "":                                       # reached a 0: the alignment starts here
            alignments.append((i, j, row_a, row_b))
        if "↖" in arrows: walk_back(i-1, j-1, a[i-1] + row_a, b[j-1] + row_b)
        if "↑" in arrows: walk_back(i-1, j,   a[i-1] + row_a, "-"    + row_b)
        if "←" in arrows: walk_back(i,   j-1, "-"    + row_a, b[j-1] + row_b)

    highest = max(max(row) for row in scores)
    if highest > 0:                                            # all 0 means nothing matches: no alignment
        for i in range(n + 1):
            for j in range(m + 1):
                if scores[i][j] == highest:
                    walk_back(i, j, "", "")
    return traceback, alignments


def show(a, b, traceback, alignments):
    print("Traceback matrix:   ↖ aligned pair   ↑ gap in b   ← gap in a   · start")
    print("    " + "".join(f"{letter:<4}" for letter in "-" + b))
    for letter, row in zip("-" + a, traceback):
        cells = [arrows or "·" for arrows in row]
        print(f"{letter}   " + "".join(f"{cell:<4}" for cell in cells))

    print(f"\n{len(alignments)} optimal alignment(s)")
    for number, (i, j, row_a, row_b) in enumerate(alignments, 1):
        used_a = row_a.replace("-", "")                        # the letters of a / b this alignment covers
        used_b = row_b.replace("-", "")
        bars = ""                                              # | marks each matching column
        for x, y in zip(row_a, row_b):
            bars += "|" if x == y else " "
        print(f"\n#{number}   a[{i}:{i + len(used_a)}]   b[{j}:{j + len(used_b)}]\n{row_a}\n{bars}\n{row_b}")


if __name__ == "__main__":
    a = input("Sequence 1: ").strip().upper()
    b = input("Sequence 2: ").strip().upper()
    match, mismatch, gap = map(int, input("Match, mismatch, gap scores (e.g. 1 -1 -2): ").split())
    traceback, alignments = smith_waterman(a, b, match, mismatch, gap)
    show(a, b, traceback, alignments)
