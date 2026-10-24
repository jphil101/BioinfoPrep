import io, os, runpy
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


def tree_sequences(edges, L=40):
    """Sequences whose p-distances equal path lengths on a known tree.
    edges = {leaves below edge: length}; each unit of length = one site where those leaves differ."""
    leaves = sorted({n for e in edges for n in e})
    cols = [set(e) for e, k in edges.items() for _ in range(k)]
    cols += [set()] * (L - len(cols))                         # pad with constant sites
    return [f"{n} {''.join('T' if n in c else 'A' for c in cols)}" for n in leaves], len(cols)
