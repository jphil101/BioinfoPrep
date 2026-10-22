import warnings; warnings.filterwarnings("ignore")   # log(0) = -inf in matrices
from Bio import Phylo
from io import StringIO
from needleman_wunsch import needleman_wunsch
from smith_waterman import smith_waterman
from distances import distances
from upgma import upgma
from neighbor_joining import neighbor_joining
from fitch_margoliash import fitch_margoliash
from blosum import blosum
from pam import pam

for align in (needleman_wunsch, smith_waterman):
    best = align("GATTACA", "GCATGCU")
    print(f"{align.__name__}  score = {best.score}\n{best}")

D = distances({"A": "ACGTACGTAC", "B": "ACGTACGTTC", "C": "ACGAACGATC", "D": "TCGAACCATC", "E": "TCGTTCCATG"})
for name, tree in [("UPGMA", upgma(D)), ("NJ", neighbor_joining(D)),
                   ("Fitch-Margoliash", Phylo.read(StringIO(fitch_margoliash(D)), "newick"))]:
    print(name); Phylo.draw_ascii(tree)

print("BLOSUM", blosum(["ACDA", "ACEA", "SCDA", "ACDG"]), sep="\n")
print("PAM250", pam(["ACDEAKLM", "ACDEAKLV", "ACEEAKLM", "SCDEAKIM"]), sep="\n")
