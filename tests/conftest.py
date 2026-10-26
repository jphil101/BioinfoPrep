import io, os, re, runpy
os.environ["MPLBACKEND"] = "Agg"                      # plots render off-screen in tests
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock
import pytest

ROOT = Path(__file__).parent.parent


@pytest.fixture
def run():
    """Run a script as if a user typed `answers` at its input() prompts.
    Returns its variables, plus what it printed under "_printed"."""
    def _run(script, *answers):
        typed, out = iter(answers), io.StringIO()
        with mock.patch("builtins.input", lambda prompt="": next(typed)), redirect_stdout(out):
            g = runpy.run_path(str(ROOT / script))
        return {**g, "_printed": out.getvalue()}
    return _run


@pytest.fixture
def run_with():
    """Run a script with its in-code `seqs = [...]` or `matrix = [[...]]` swapped for test data.
    Returns its variables, plus what it printed under "_printed"."""
    def _run(script, seqs=None, matrix=None):
        import matplotlib.pyplot as plt
        plt.close("all")                                     # keep only this run's plot open
        src = (ROOT / script).read_text()
        if seqs is not None:
            src, n = re.subn(r"^seqs = \[.*?\]", f"seqs = {seqs!r}", src, count=1, flags=re.M)
            assert n == 1, "no seqs = [...] line found"
        if matrix is not None:
            src, n = re.subn(r"^matrix = \[\[.*?\]\]", f"matrix = {matrix!r}", src, count=1, flags=re.M | re.S)
            assert n == 1, "no matrix = [[...]] block found"
        g, out = {"__name__": "__main__"}, io.StringIO()
        with redirect_stdout(out):
            exec(compile(src, str(ROOT / script), "exec"), g)
        return {**g, "_printed": out.getvalue()}
    return _run


def tree_sequences(edges, L=40):
    """Sequences whose p-distances equal path lengths on a known tree.
    edges = {leaves below edge: length}; each unit of length = one site where those leaves differ."""
    leaves = sorted({n for e in edges for n in e})
    cols = [set(e) for e, k in edges.items() for _ in range(k)]
    cols += [set()] * (L - len(cols))                         # pad with constant sites
    return [''.join('T' if n in c else 'A' for c in cols) for n in leaves], leaves
