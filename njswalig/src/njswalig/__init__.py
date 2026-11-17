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

Built the score matrix yourself? ``njalign.traceback(seq1, seq2, matrix, scores)`` (or
``swalign.traceback``) checks it, and returns the traceback matrix and every optimal alignment.
"""

from itertools import combinations

__all__ = ["njalign", "swalign", "Aligner", "Alignments", "PairAlignment", "Traceback"]
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
        header = f"{_MODES[self.mode]}: {self.name1} vs {self.name2}, score = {_num(self.score)}"
        if not self.aligned1:
            return header + "\n(no positive-scoring local alignment)"
        return "\n".join([header] + self._rows())

    def _rows(self):
        """The three printed rows: sequence 1, the match line, sequence 2."""
        w = max(len(self.name1), len(self.name2))
        if self.mode == "global":
            return [f"{self.name1:<{w}}  {self.aligned1}",
                    f"{'':<{w}}  {self.match_line}",
                    f"{self.name2:<{w}}  {self.aligned2}"]
        p = len(str(max(self.start1, self.start2)))    # local: show where the region sits in each sequence
        return [f"{self.name1:<{w}}  {self.start1:>{p}} {self.aligned1} {self.end1}",
                f"{'':<{w}}  {'':>{p}} {self.match_line}",
                f"{self.name2:<{w}}  {self.start2:>{p}} {self.aligned2} {self.end2}"]

    def matrix_str(self):
        """The DP score table with sequence letters along the edges."""
        return _table(self.seq1, self.seq2, [[_num(v) for v in row] for row in self.matrix])

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


class Traceback:
    """What :meth:`Aligner.traceback` returns: everything you can read off a filled score matrix.

    Attributes:

    - ``matrix``: the score matrix you passed in (rows follow seq1, columns follow seq2).
    - ``arrows``: the traceback matrix, same shape, as strings made of ``↖`` (seq1 letter paired
      with seq2 letter), ``↑`` (seq1 letter against a gap) and ``←`` (gap against a seq2 letter);
      ``""`` where no arrow leaves the cell. Ties give several arrows, in the order ``↖↑←``.
    - ``alignments``: every optimal alignment, as :class:`PairAlignment` objects (``score``,
      ``aligned1``, ``aligned2``, ``match_line``, ``start1``/``end1``, ``start2``/``end2``, ...).
    - ``score``: the optimal score. ``complete``: False if ``limit`` stopped the search early.
    - ``mode``, ``seq1``, ``seq2``, ``name1``, ``name2``, ``scores`` (``(match, mismatch, gap)``).
    """

    def __init__(self, mode, name1, name2, seq1, seq2, scores, matrix, arrows, alignments, complete):
        self.mode, self.name1, self.name2, self.seq1, self.seq2 = mode, name1, name2, seq1, seq2
        self.scores, self.matrix, self.arrows = scores, matrix, arrows
        self.alignments, self.complete = alignments, complete
        if mode == "global":
            self.score = matrix[-1][-1]
        else:
            self.score = max(max(row) for row in matrix)

    def __str__(self):
        match, mismatch, gap = self.scores
        lines = [f"{_MODES[self.mode]}: {self.name1} vs {self.name2}, "
                 f"match {_num(match)}, mismatch {_num(mismatch)}, gap {_num(gap)}",
                 "",
                 "Score matrix:",
                 _table(self.seq1, self.seq2, [[_num(v) for v in row] for row in self.matrix]),
                 "",
                 "Traceback matrix (↖ letter with letter, ↑ seq1 letter with a gap, ← seq2 letter with a gap):",
                 _table(self.seq1, self.seq2, [[cell or "·" for cell in row] for row in self.arrows]),
                 ""]
        count = len(self.alignments)
        if count == 0:
            lines.append("No local alignment: no cell scores above 0.")
            return "\n".join(lines)
        lines.append(f"{count} optimal alignment{'' if count == 1 else 's'}, score {_num(self.score)}:")
        for number, pair in enumerate(self.alignments, 1):
            lines += ["", f"#{number}  score = {_num(pair.score)}"] + pair._rows()
        if not self.complete:
            lines += ["", f"(stopped after {count} alignments because of limit={count}; there are more)"]
        return "\n".join(lines)

    def show(self):
        """Print the score matrix, the traceback matrix, then every optimal alignment with its score."""
        print(self)

    def __repr__(self):
        return (f"<Traceback {self.mode} {self.name1} vs {self.name2} score={_num(self.score)} "
                f"alignments={len(self.alignments)}{'' if self.complete else '+'}>")


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

    def traceback(self, seq1, seq2, matrix, scores, limit=None, names=None):
        """Trace back through a score matrix you filled in yourself.

        ``matrix[i][j]`` is the score for the first ``i`` letters of ``seq1`` against the first
        ``j`` letters of ``seq2``, so it has ``len(seq1)+1`` rows of ``len(seq2)+1`` numbers,
        border row and column included. Every cell is checked against ``scores = [match,
        mismatch, gap]``; the first wrong one raises ``ValueError`` naming the cell.

        Returns a :class:`Traceback` holding the traceback matrix (``.arrows``) and every optimal
        alignment (``.alignments``); ``.show()`` prints them. The number of tied alignments can
        grow very fast (e.g. with gap 0); ``limit`` stops after that many.
        """
        if not isinstance(seq1, str) or not isinstance(seq2, str):
            raise ValueError("seq1 and seq2 must be strings, e.g. 'ATCG'")
        (a, b), (name1, name2) = _check_seqs([seq1, seq2], names)
        match, mismatch, gap = _check_scores(scores)
        if limit is not None and (not isinstance(limit, int) or isinstance(limit, bool) or limit < 1):
            raise ValueError(f"limit must be a whole number of at least 1, or None; got {limit!r}")
        S = _check_matrix(self.mode, a, b, matrix, match, mismatch, gap)
        arrows = _arrows(self.mode, a, b, S, match, mismatch, gap)
        found, complete = _walk(self.mode, a, b, S, arrows, limit)
        alignments = [PairAlignment(self.mode, name1, name2, a, b, row1, row2,
                                    _column_score(row1, row2, match, mismatch, gap), S, start, end)
                      for start, end, row1, row2 in found]
        return Traceback(self.mode, name1, name2, a, b, (match, mismatch, gap), S, arrows, alignments, complete)

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


def _check_matrix(mode, a, b, matrix, match, mismatch, gap):
    """Check the shape and every cell of a hand-built score matrix; return a copy as a list of lists."""
    rows, cols = len(a) + 1, len(b) + 1
    if not isinstance(matrix, (list, tuple)) or not all(isinstance(row, (list, tuple)) for row in matrix):
        raise ValueError("matrix must be a list of rows, each a list of numbers")
    for row in matrix:
        for value in row:
            if not isinstance(value, (int, float)) or isinstance(value, bool):
                raise ValueError(f"matrix must hold only numbers; found {value!r}")
    if len(matrix) != rows or any(len(row) != cols for row in matrix):
        got = f"{len(matrix)} rows of {sorted({len(row) for row in matrix}) or [0]} numbers"
        if len(matrix) == cols and all(len(row) == rows for row in matrix):
            raise ValueError(f"matrix is transposed: it needs {rows} rows (one per letter of seq1 "
                             f"{a!r}, plus the border row) of {cols} numbers (one per letter of seq2 "
                             f"{b!r}, plus the border column), but has {cols} rows of {rows}. "
                             f"Use matrix[i][j] with i for seq1 and j for seq2.")
        raise ValueError(f"matrix must have {rows} rows (len(seq1) + 1, seq1 = {a!r}) of {cols} numbers "
                         f"(len(seq2) + 1, seq2 = {b!r}); got {got}")

    local = mode == "local"
    for i in range(rows):
        for j in range(cols):
            if i == 0 or j == 0:                       # border: global pays a gap per letter, local starts free
                expected = 0 if local else (i + j) * gap
            else:
                pair = match if a[i - 1] == b[j - 1] else mismatch
                expected = max(matrix[i - 1][j - 1] + pair, matrix[i - 1][j] + gap, matrix[i][j - 1] + gap)
                if local:
                    expected = max(expected, 0)
            if not _same(matrix[i][j], expected):
                where = []
                if i:
                    where.append(f"row {i} is seq1 letter {i} {a[i - 1]!r}")
                if j:
                    where.append(f"column {j} is seq2 letter {j} {b[j - 1]!r}")
                raise ValueError(f"cell ({i}, {j}) is {_num(matrix[i][j])} but the scores give {_num(expected)}"
                                 + (f" ({', '.join(where)})" if where else ""))
    return [list(row) for row in matrix]


def _arrows(mode, a, b, S, match, mismatch, gap):
    """The traceback matrix: from each cell, every step back that explains its score."""
    arrows = []
    for i in range(len(a) + 1):
        row = []
        for j in range(len(b) + 1):
            cell = ""
            if not (mode == "local" and _same(S[i][j], 0)):   # local: a 0 cell starts an alignment, no arrow
                if i > 0 and j > 0 and _same(S[i][j], S[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch)):
                    cell += "↖"
                if i > 0 and _same(S[i][j], S[i - 1][j] + gap):
                    cell += "↑"
                if j > 0 and _same(S[i][j], S[i][j - 1] + gap):
                    cell += "←"
            row.append(cell)
        arrows.append(row)
    return arrows


def _walk(mode, a, b, S, arrows, limit):
    """Follow every arrow back from the end cell(s). Returns ([(start, end, row1, row2), ...], complete).

    Uses an explicit stack (no recursion), so long sequences are fine. The rows are built from the
    end of the alignment, so each new letter goes on the front.
    """
    local = mode == "local"
    if local:
        best = max(max(row) for row in S)
        if best <= 0 or _same(best, 0):                # nothing scores above 0: no local alignment
            return [], True
        ends = [(i, j) for i in range(len(a) + 1) for j in range(len(b) + 1) if _same(S[i][j], best)]
    else:
        ends = [(len(a), len(b))]

    found = []
    for number, end in enumerate(ends):
        stack = [(end[0], end[1], "", "")]
        while stack:
            i, j, row1, row2 = stack.pop()
            if (local and _same(S[i][j], 0)) or (not local and i == 0 and j == 0):
                found.append(((i, j), end, row1, row2))     # reached the start of this alignment
                if limit is not None and len(found) >= limit:
                    return found, not stack and number == len(ends) - 1
                continue                                    # never step past the start
            for arrow in reversed(arrows[i][j]):            # reversed, so ↖ is explored first
                if arrow == "↖":
                    stack.append((i - 1, j - 1, a[i - 1] + row1, b[j - 1] + row2))
                elif arrow == "↑":
                    stack.append((i - 1, j, a[i - 1] + row1, "-" + row2))
                else:
                    stack.append((i, j - 1, "-" + row1, b[j - 1] + row2))
    return found, True


def _column_score(row1, row2, match, mismatch, gap):
    """Total score of an alignment, column by column."""
    total = 0
    for x, y in zip(row1, row2):
        total += gap if "-" in (x, y) else match if x == y else mismatch
    return total


def _table(seq1, seq2, cells):
    """Cells (strings) laid out under seq2's letters, with seq1's letters down the left."""
    w = max(len(c) for row in cells for c in row)
    lines = ["  " + " ".join(f"{c:>{w}}" for c in "-" + seq2)]
    for letter, row in zip("-" + seq1, cells):
        lines.append(f"{letter} " + " ".join(f"{c:>{w}}" for c in row))
    return "\n".join(lines)


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
