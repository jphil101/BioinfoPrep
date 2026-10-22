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


# join closest pair, new node at height d/2, repeat
# (Biopython's .upgma() is really WPGMA, so it's written out here)
def upgma(D):
    d = lambda X, Y: np.mean([D[i, j] for i in X for j in Y])     # avg leaf-to-leaf distance
    C = {(n,): (n, 0) for n in D.names}                            # cluster -> (newick, height)
    while len(C) > 1:
        A, B = min(combinations(C, 2), key=lambda p: d(*p))
        h = d(A, B) / 2
        (na, ha), (nb, hb) = C.pop(A), C.pop(B)
        C[A + B] = (f"({na}:{h-ha},{nb}:{h-hb})", h)
    (newick, _), = C.values()
    return Phylo.read(StringIO(newick + ";"), "newick")


if __name__ == "__main__":
    seqs = {"A": "ACGTACGTAC", "B": "ACGTACGTTC", "C": "ACGAACGATC", "D": "TCGAACCATC", "E": "TCGTTCCATG"}
    Phylo.draw_ascii(upgma(distances(seqs)))
