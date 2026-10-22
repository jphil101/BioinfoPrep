from Bio.Phylo.TreeConstruction import DistanceMatrix

def distances(s):             # {name: aligned seq} -> p-distance (fraction of mismatches)
    n = list(s)
    return DistanceMatrix(n, [[sum(x != y for x, y in zip(s[a], s[b])) / len(s[a]) for b in n[:i+1]] for i, a in enumerate(n)])
