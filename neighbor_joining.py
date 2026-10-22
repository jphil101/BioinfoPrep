from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# {name: aligned seq} -> p-distance matrix (1 - identity)
def distances(seqs):
    aln = MultipleSeqAlignment([SeqRecord(Seq(s), id=n) for n, s in seqs.items()])
    return DistanceCalculator("identity").get_distance(aln)


def neighbor_joining(D):
    return DistanceTreeConstructor().nj(D)


if __name__ == "__main__":
    seqs = {"A": "ACGTACGTAC", "B": "ACGTACGTTC", "C": "ACGAACGATC", "D": "TCGAACCATC", "E": "TCGTTCCATG"}
    Phylo.draw_ascii(neighbor_joining(distances(seqs)))
