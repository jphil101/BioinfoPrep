"""Needleman-Wunsch (global) and Smith-Waterman (local) alignment. No imports, linear gap penalty.

    python seqal.py                                  # interactive
    trace, alignments = align("GATTACA", "GCATGCU", local=False, match=1, mismatch=-1, gap=-1)
"""

STEP = {"D": (1, 1), "U": (1, 0), "L": (0, 1)}   # arrow -> (steps in a, steps in b)


def align(a, b, local=False, match=1, mismatch=-1, gap=-2):
    """local=False: Needleman-Wunsch.  local=True: Smith-Waterman.  Use integer scores (exact ties).

    Returns (trace, alignments):
      trace[i][j]  arrows into cell (i, j), one letter per tied route:
                   D diagonal (a[i-1] with b[j-1]), U up (a[i-1] with a gap), L left (gap with b[j-1]).
                   "" means an alignment starts here.
      alignments   every optimal alignment as (start_i, start_j, row_a, row_b); "-" is a gap.
    """
    n, m = len(a), len(b)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[""] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1):
        for j in range(m + 1):
            moves = {}                                    # arrow -> score we would get by taking it
            if i and j: moves["D"] = score[i-1][j-1] + (match if a[i-1] == b[j-1] else mismatch)
            if i:       moves["U"] = score[i-1][j] + gap
            if j:       moves["L"] = score[i][j-1] + gap
            if not moves: continue                        # cell (0, 0)
            best = max(moves.values())
            if local and best <= 0: continue              # SW: scores never drop below 0, alignment starts here
            score[i][j] = best
            trace[i][j] = "".join(k for k in moves if moves[k] == best)

    def trace_back(i, j):
        """Yield every optimal alignment that ends at (i, j), following every tied arrow."""
        if not trace[i][j]:
            yield i, j, "", ""                            # reached the start
        for k in trace[i][j]:
            di, dj = STEP[k]
            for si, sj, x, y in trace_back(i - di, j - dj):
                yield si, sj, x + (a[i-1] if di else "-"), y + (b[j-1] if dj else "-")

    if local:                                             # SW: start from every cell holding the best score
        top = max(map(max, score))
        ends = [(i, j) for i in range(n + 1) for j in range(m + 1) if top > 0 and score[i][j] == top]
    else:                                                 # NW: start from the bottom-right corner
        ends = [(n, m)]
    return trace, [al for i, j in ends for al in trace_back(i, j)]


def show(a, b, trace, alignments):
    arrows = {"D": "↖", "U": "↑", "L": "←"}
    cells = [["".join(arrows[k] for k in t) or "·" for t in row] for row in trace]
    w = max(len(c) for row in cells for c in row)
    print("Traceback matrix:")
    print(" ", *(f"{c:>{w}}" for c in " " + b))
    for ch, row in zip(" " + a, cells):
        print(ch, *(f"{c:>{w}}" for c in row))
    print(f"\n{len(alignments)} optimal alignment(s)")
    for k, (i, j, x, y) in enumerate(alignments, 1):
        bars = "".join("|" if p == q else " " for p, q in zip(x, y))
        print(f"\n#{k}  a[{i}:{i + len(x) - x.count('-')}]  b[{j}:{j + len(y) - y.count('-')}]\n{x}\n{bars}\n{y}")


if __name__ == "__main__":
    a = input("Sequence 1: ").strip().upper()
    b = input("Sequence 2: ").strip().upper()
    local = input("Global (NW) or local (SW)? [g/l]: ").strip().lower() == "l"
    match, mismatch, gap = map(int, input("Match, mismatch, gap scores (e.g. 1 -1 -2): ").split())
    show(a, b, *align(a, b, local, match, mismatch, gap))
