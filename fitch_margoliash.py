import numpy as np
from io import StringIO
from itertools import combinations
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# {name: aligned seq} -> p-distance matrix (1 - identity)
def distances(seqs):
    aln = MultipleSeqAlignment([SeqRecord(Seq(s), id=n) for n, s in seqs.items()])
    return DistanceCalculator("identity").get_distance(aln)


# join closest pair, 3-point formula vs. the rest, repeat
def fitch_margoliash(D):
    d = lambda X, Y: np.mean([D[i, j] for i in X for j in Y])     # avg leaf-to-leaf distance
    C = {(n,): (n, 0) for n in D.names}                            # cluster -> (newick, avg depth)
    while len(C) > 2:
        A, B = min(combinations(C, 2), key=lambda p: d(*p))
        R = sum((c for c in C if c not in (A, B)), ())             # all other leaves
        a = (d(A, B) + d(A, R) - d(B, R)) / 2
        b = d(A, B) - a
        (na, ha), (nb, hb) = C.pop(A), C.pop(B)
        C[A + B] = (f"({na}:{a-ha},{nb}:{b-hb})", (len(A)*a + len(B)*b) / len(A+B))
    (X, (nx, hx)), (Y, (ny, hy)) = C.items()
    return Phylo.read(StringIO(f"({nx},{ny}:{d(X, Y)-hx-hy});"), "newick")


if __name__ == "__main__":
    seqs = {"A": "ACGTACGTAC", "B": "ACGTACGTTC", "C": "ACGAACGATC", "D": "TCGAACCATC", "E": "TCGTTCCATG"}
    Phylo.draw_ascii(fitch_margoliash(distances(seqs)))
