"""Needleman-Wunsch (global alignment). No imports, linear gap penalty, integer scores."""

STEP = {"↖": (1, 1), "↑": (1, 0), "←": (0, 1)}      # arrow -> (steps in a, steps in b)


def needleman_wunsch(a, b, match=1, mismatch=-1, gap=-2):
    """Returns (trace, alignments).
    trace[i][j]  every tied arrow into cell (i, j): ↖ a[i-1] with b[j-1], ↑ a[i-1] with a gap, ← gap with b[j-1]
    alignments   every optimal alignment as (row_a, row_b); '-' is a gap"""
    n, m = len(a), len(b)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[""] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1): score[i][0], trace[i][0] = i * gap, "↑"     # a against nothing but gaps
    for j in range(1, m + 1): score[0][j], trace[0][j] = j * gap, "←"     # b against nothing but gaps

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            moves = {"↖": score[i-1][j-1] + (match if a[i-1] == b[j-1] else mismatch),
                     "↑": score[i-1][j] + gap,
                     "←": score[i][j-1] + gap}
            score[i][j] = max(moves.values())
            trace[i][j] = "".join(k for k in moves if moves[k] == score[i][j])     # keep every tie

    def paths(i, j):                                  # every optimal alignment of a[:i] with b[:j]
        if not trace[i][j]: yield "", ""              # reached cell (0, 0)
        for k in trace[i][j]:
            di, dj = STEP[k]
            for x, y in paths(i - di, j - dj):
                yield x + (a[i-1] if di else "-"), y + (b[j-1] if dj else "-")

    return trace, list(paths(n, m))                   # global: trace back from the bottom-right corner


def show(a, b, trace, alignments):
    print("Traceback matrix:")
    print(" ", *(f"{c:<3}" for c in " " + b))
    for ch, row in zip(" " + a, trace):
        print(ch, *(f"{t or '·':<3}" for t in row))
    print(f"\n{len(alignments)} optimal alignment(s)")
    for k, (x, y) in enumerate(alignments, 1):
        print(f"\n#{k}\n{x}\n{''.join('|' if p == q else ' ' for p, q in zip(x, y))}\n{y}")


if __name__ == "__main__":
    a, b = input("Sequence 1: ").strip().upper(), input("Sequence 2: ").strip().upper()
    scores = map(int, input("Match, mismatch, gap scores (e.g. 1 -1 -2): ").split())
    show(a, b, *needleman_wunsch(a, b, *scores))
