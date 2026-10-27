"""njswalig -- Needleman-Wunsch (global) and Smith-Waterman (local) pairwise alignment.

    >>> from njswalig import njalign, swalign
    >>> seqs = ["ATCG", "TTTA"]
    >>> scores = [2, -1, -2]                 # [match, mismatch, gap]
    >>> alignment = njalign.sequences(seqs, scores)
    >>> alignment.show()
    Needleman-Wunsch (global): Seq1 vs Seq2, score = -1
    Seq1  ATCG
          .|..
    Seq2  TTTA

Every pair of sequences in the list is aligned. ``sequence`` is the same as ``sequences``.
Scoring is linear: each gap position costs ``gap``.
"""

from itertools import combinations

__all__ = ["njalign", "swalign", "Aligner", "Alignments", "PairAlignment"]
__version__ = "0.1.0"

_MODES = {"global": "Needleman-Wunsch (global)", "local": "Smith-Waterman (local)"}


class PairAlignment:
    """The best alignment of one pair of sequences.

    Attributes: ``name1``, ``name2``, ``seq1``, ``seq2``, ``aligned1``, ``aligned2`` (with ``-``
    for gaps), ``score``, ``mode`` ("global"/"local"), ``matrix`` (the filled DP score table,
    ``len(seq1)+1`` rows by ``len(seq2)+1`` columns), and 1-based ``start1``/``end1``,
    ``start2``/``end2`` giving the aligned region of each sequence.
    """

    def __init__(self, mode, name1, name2, seq1, seq2, aligned1, aligned2, score, matrix, start, end):
        self.mode, self.name1, self.name2, self.seq1, self.seq2 = mode, name1, name2, seq1, seq2
        self.aligned1, self.aligned2, self.score, self.matrix = aligned1, aligned2, score, matrix
        self.start1, self.start2 = start[0] + 1, start[1] + 1
        self.end1, self.end2 = end

    @property
    def match_line(self):
        """``|`` match, ``.`` mismatch, space for a gap."""
        return "".join(" " if "-" in (x, y) else "|" if x == y else "." for x, y in zip(self.aligned1, self.aligned2))

    def __str__(self):
        lines = [f"{_MODES[self.mode]}: {self.name1} vs {self.name2}, score = {_num(self.score)}"]
        if not self.aligned1:
            return lines[0] + "\n(no positive-scoring local alignment)"
        w = max(len(self.name1), len(self.name2))
        if self.mode == "global":
            lines += [f"{self.name1:<{w}}  {self.aligned1}",
                      f"{'':<{w}}  {self.match_line}",
                      f"{self.name2:<{w}}  {self.aligned2}"]
        else:                                          # local: show where the region sits in each sequence
            p = len(str(max(self.start1, self.start2)))
            lines += [f"{self.name1:<{w}}  {self.start1:>{p}} {self.aligned1} {self.end1}",
                      f"{'':<{w}}  {'':>{p}} {self.match_line}",
                      f"{self.name2:<{w}}  {self.start2:>{p}} {self.aligned2} {self.end2}"]
        return "\n".join(lines)

    def matrix_str(self):
        """The DP score table with sequence letters along the edges."""
        cells = [[_num(v) for v in row] for row in self.matrix]
        w = max(len(c) for row in cells for c in row)
        head = "  " + " ".join(f"{c:>{w}}" for c in "-" + self.seq2)
        return "\n".join([head] + [f"{r} " + " ".join(f"{c:>{w}}" for c in row)
                                   for r, row in zip("-" + self.seq1, cells)])

    def show(self, matrix=False):
        """Print the alignment (and the DP score table if ``matrix=True``)."""
        print(self)
        if matrix:
            print(self.matrix_str())

    def __repr__(self):
        return f"<PairAlignment {self.mode} {self.name1} vs {self.name2} score={_num(self.score)}>"


class Alignments(list):
    """A list of :class:`PairAlignment`, one per pair of input sequences."""

    def show(self, matrix=False):
        """Print every alignment (and DP score tables if ``matrix=True``), separated by blank lines."""
        for i, pair in enumerate(self):
            if i:
                print()
            pair.show(matrix=matrix)

    def __str__(self):
        return "\n\n".join(map(str, self))


class Aligner:
    """``Aligner("global")`` is Needleman-Wunsch; ``Aligner("local")`` is Smith-Waterman."""

    def __init__(self, mode):
        if mode not in _MODES:
            raise ValueError(f"mode must be 'global' or 'local', got {mode!r}")
        self.mode = mode

    def sequences(self, seqs, scores, names=None):
        """Align every pair in ``seqs`` using ``scores = [match, mismatch, gap]``.

        Sequences are upper-cased and stripped of surrounding whitespace. Names default to
        Seq1, Seq2, ... in list order. Returns :class:`Alignments` (call ``.show()`` to print).
        """
        seqs, names = _check_seqs(seqs, names)
        match, mismatch, gap = _check_scores(scores)
        return Alignments(_align(self.mode, names[i], names[j], seqs[i], seqs[j], match, mismatch, gap)
                          for i, j in combinations(range(len(seqs)), 2))

    sequence = sequences                               # same thing, either name works

    def __repr__(self):
        return f"<Aligner {_MODES[self.mode]}>"


njalign = Aligner("global")
swalign = Aligner("local")


def _align(mode, name1, name2, a, b, match, mismatch, gap):
    local = mode == "local"
    n, m = len(a), len(b)
    S = [[0] * (m + 1) for _ in range(n + 1)]
    if not local:                                      # global: leading gaps cost gap each
        for i in range(1, n + 1):
            S[i][0] = i * gap
        for j in range(1, m + 1):
            S[0][j] = j * gap

    def pair_score(i, j):
        return match if a[i - 1] == b[j - 1] else mismatch

    best = (0, 0, 0)                                   # local: (score, i, j) of the highest cell
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            v = max(S[i - 1][j - 1] + pair_score(i, j), S[i - 1][j] + gap, S[i][j - 1] + gap)
            if local:
                v = max(v, 0)
                if v > best[0]:
                    best = (v, i, j)
            S[i][j] = v

    i, j = (best[1], best[2]) if local else (n, m)
    end, score = (i, j), (best[0] if local else S[n][m])
    out1, out2 = [], []
    while i > 0 or j > 0:                              # traceback, preferring diagonal, then up, then left
        if local and S[i][j] == 0:
            break
        if i > 0 and j > 0 and _same(S[i][j], S[i - 1][j - 1] + pair_score(i, j)):
            out1.append(a[i - 1]); out2.append(b[j - 1]); i -= 1; j -= 1
        elif i > 0 and _same(S[i][j], S[i - 1][j] + gap):
            out1.append(a[i - 1]); out2.append("-"); i -= 1
        else:
            out1.append("-"); out2.append(b[j - 1]); j -= 1
    return PairAlignment(mode, name1, name2, a, b, "".join(reversed(out1)), "".join(reversed(out2)),
                         score, S, (i, j), end)


def _same(x, y):
    return abs(x - y) < 1e-9


def _num(x):
    return str(int(x)) if float(x).is_integer() else f"{x:g}"


def _check_scores(scores):
    ok = (isinstance(scores, (list, tuple)) and len(scores) == 3
          and all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in scores))
    if not ok:
        raise ValueError(f"scores must be [match, mismatch, gap], e.g. [2, -1, -2]; got {scores!r}")
    return scores


def _check_seqs(seqs, names):
    if isinstance(seqs, str) or not isinstance(seqs, (list, tuple)) or not all(isinstance(s, str) for s in seqs):
        raise ValueError("seqs must be a list of sequences (strings), e.g. ['ATCG', 'TTTA']")
    if len(seqs) < 2:
        raise ValueError("need at least 2 sequences to align")
    seqs = [s.strip().upper() for s in seqs]
    names = list(names) if names is not None else [f"Seq{i + 1}" for i in range(len(seqs))]
    if len(names) != len(seqs) or len(set(names)) != len(names):
        raise ValueError(f"need {len(seqs)} unique names, got {names!r}")
    return seqs, names
