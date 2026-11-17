"""Smith-Waterman (local alignment). No imports, linear gap penalty, integer scores."""

STEP = {"↖": (1, 1), "↑": (1, 0), "←": (0, 1)}      # arrow -> (steps in a, steps in b)


def smith_waterman(a, b, match=1, mismatch=-1, gap=-2):
    """Returns (trace, alignments).
    trace[i][j]  every tied arrow into cell (i, j): ↖ a[i-1] with b[j-1], ↑ a[i-1] with a gap, ← gap with b[j-1];
                 "" means an alignment starts here
    alignments   every optimal alignment as (start_i, start_j, row_a, row_b); '-' is a gap"""
    n, m = len(a), len(b)
    score = [[0] * (m + 1) for _ in range(n + 1)]     # borders stay 0: no gap penalty before the alignment starts
    trace = [[""] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            moves = {"↖": score[i-1][j-1] + (match if a[i-1] == b[j-1] else mismatch),
                     "↑": score[i-1][j] + gap,
                     "←": score[i][j-1] + gap}
            best = max(moves.values())
            if best > 0:                              # otherwise the score stays 0 and an alignment can start here
                score[i][j] = best
                trace[i][j] = "".join(k for k in moves if moves[k] == best)            # keep every tie

    def paths(i, j):                                  # every optimal alignment ending at (i, j)
        if not trace[i][j]: yield i, j, "", ""        # reached a 0: this is where it starts
        for k in trace[i][j]:
            di, dj = STEP[k]
            for si, sj, x, y in paths(i - di, j - dj):
                yield si, sj, x + (a[i-1] if di else "-"), y + (b[j-1] if dj else "-")

    top = max(map(max, score))                        # local: trace back from every cell holding the best score
    ends = [(i, j) for i in range(n + 1) for j in range(m + 1) if top > 0 and score[i][j] == top]
    return trace, [al for i, j in ends for al in paths(i, j)]


def show(a, b, trace, alignments):
    print("Traceback matrix:")
    print(" ", *(f"{c:<3}" for c in " " + b))
    for ch, row in zip(" " + a, trace):
        print(ch, *(f"{t or '·':<3}" for t in row))
    print(f"\n{len(alignments)} optimal alignment(s)")
    for k, (i, j, x, y) in enumerate(alignments, 1):
        bars = "".join("|" if p == q else " " for p, q in zip(x, y))
        print(f"\n#{k}  a[{i}:{i + len(x) - x.count('-')}]  b[{j}:{j + len(y) - y.count('-')}]\n{x}\n{bars}\n{y}")


if __name__ == "__main__":
    a, b = input("Sequence 1: ").strip().upper(), input("Sequence 2: ").strip().upper()
    scores = map(int, input("Match, mismatch, gap scores (e.g. 1 -1 -2): ").split())
    show(a, b, *smith_waterman(a, b, *scores))
