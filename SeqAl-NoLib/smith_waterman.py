"""Smith-Waterman: LOCAL alignment (finds the best-matching region; the rest of each sequence is ignored).

No imports. Scores are whole numbers; `gap` is the cost of one gap position, so use a negative number.

    python smith_waterman.py
"""


def smith_waterman(a, b, match=1, mismatch=-1, gap=-2):
    """Returns (scores, alignments).

    scores[i][j] is the best score of an alignment ENDING at a[i-1] and b[j-1] (never below 0).
    alignments is a list of every optimal alignment as (start_i, start_j, row_a, row_b).
    It begins at a[start_i] and b[start_j]; "-" marks a gap.
    """
    n, m = len(a), len(b)

    def pair_score(i, j):                                      # lining up a[i-1] with b[j-1]
        return match if a[i-1] == b[j-1] else mismatch

    # 1. Table with one extra row and column for "nothing aligned yet".
    #    Unlike Needleman-Wunsch, the first row and column stay 0: skipping the start of a sequence is free.
    scores = [[0] * (m + 1) for _ in range(n + 1)]

    # 2. Every other cell: the best of the three ways to arrive there, but never below 0.
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = scores[i-1][j-1] + pair_score(i, j)     # a[i-1] with b[j-1]
            up = scores[i-1][j] + gap                          # a[i-1] with a gap
            left = scores[i][j-1] + gap                        # b[j-1] with a gap
            scores[i][j] = max(diagonal, up, left, 0)          # 0 means "start a fresh alignment here"

    # 3. Walk back from every cell holding the highest score, until we reach a 0.
    #    A step back is valid if it explains this cell's score: the neighbour's score plus that move's
    #    cost equals it. If several steps are valid (a tie), follow them all.
    #    We meet the alignment from its end, so each new letter goes on the FRONT of the rows.
    alignments = []

    def walk_back(i, j, row_a, row_b):
        if scores[i][j] == 0:                                  # reached a 0: the alignment starts here
            alignments.append((i, j, row_a, row_b))
            return
        if scores[i][j] == scores[i-1][j-1] + pair_score(i, j):
            walk_back(i-1, j-1, a[i-1] + row_a, b[j-1] + row_b)
        if scores[i][j] == scores[i-1][j] + gap:
            walk_back(i-1, j,   a[i-1] + row_a, "-"    + row_b)
        if scores[i][j] == scores[i][j-1] + gap:
            walk_back(i,   j-1, "-"    + row_a, b[j-1] + row_b)

    highest = max(max(row) for row in scores)
    if highest > 0:                                            # all 0 means nothing matches: no alignment
        for i in range(n + 1):
            for j in range(m + 1):
                if scores[i][j] == highest:
                    walk_back(i, j, "", "")
    return scores, alignments


def show(a, b, scores, alignments):
    width = max(len(str(value)) for row in scores for value in row) + 1
    print("Score matrix:")
    print("  " + "".join(f"{letter:>{width}}" for letter in "-" + b))
    for letter, row in zip("-" + a, scores):
        print(f"{letter} " + "".join(f"{value:>{width}}" for value in row))

    score = max(max(row) for row in scores)                    # every optimal alignment has this score
    print(f"\n{len(alignments)} optimal alignment(s), score {score} each")
    for number, (i, j, row_a, row_b) in enumerate(alignments, 1):
        used_a = row_a.replace("-", "")                        # the letters of a / b this alignment covers
        used_b = row_b.replace("-", "")
        bars = ""                                              # | marks each matching column
        for x, y in zip(row_a, row_b):
            bars += "|" if x == y else " "
        print(f"\n#{number}   score {score}   a[{i}:{i + len(used_a)}]   b[{j}:{j + len(used_b)}]\n{row_a}\n{bars}\n{row_b}")


if __name__ == "__main__":
    a = input("Sequence 1: ").strip().upper()
    b = input("Sequence 2: ").strip().upper()
    match, mismatch, gap = map(int, input("Match, mismatch, gap scores (e.g. 1 -1 -2): ").split())
    scores, alignments = smith_waterman(a, b, match, mismatch, gap)
    show(a, b, scores, alignments)
