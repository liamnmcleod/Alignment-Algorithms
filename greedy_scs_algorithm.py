import itertools

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

"""
Function determines the length of the longest overlap between a and b
"""
def overlap(a, b, min_length=3):

    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffix in a
        if start == -1:  # no more occurrences to right
            return 0

        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move past previous match

"""
Function finds the set of shortest common superstrings of given strings.
*Given strings must have the same length*
"""
def scs(ss):
    shortest_sup = []
    for perm in itertools.permutations(ss):
        sup = perm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(perm[i], perm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += perm[i+1][olen:]
        if len(shortest_sup) == 0 or len(sup) < len(shortest_sup[0]):
            shortest_sup = [sup]  # found shorter superstring
        elif len(sup) == len(shortest_sup[0]):
            shortest_sup.append(sup)
    return shortest_sup  # return shortest

"""
Given a set of reads, this function finds the pair which overlap the most and calculates the length of the overlap.
"""
def pick_maximal_overlap(reads, k):
    reada, readb = None, None
    best_olen = 0
    for a,b in itertools.permutations(reads, 2):
        olen = overlap(a, b, k)
        if olen > best_olen:
            reada, readb = a, b
            best_olen = olen
    return reada, readb, best_olen

"""
This function implements the greedy shortest common superstring algorithm.
"""
def greedy_scs(reads, k):
    read_a, read_b, olen = pick_maximal_overlap(reads, k)
    while olen > 0:
        print(len(reads))
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = pick_maximal_overlap(reads, k)
    return ''.join(reads)

"""
This is an accelerated version of pick_maximal_overlap(reads, k).
This is achieved by building an k-mer index so that not every permutation of reads is considered.
"""
def pick_maximal_overlap_index(reads, k):
    index = {}
    for read in reads:
        kmers = []
        for i in range(len(read) - k + 1):
            kmers.append(read[i:i+k])
        for kmer in kmers:
            if kmer not in index:
                index[kmer] = set()
            index[kmer].add(read)
    for read in reads:
        for i in range(len(read)-k+1):
            dummy = read[i:i+k]
            if dummy not in index:
                index[dummy] = set()
            index[dummy].add(read)
    reada, readb = None, None
    best_olen = 0
    for a in reads:
        for b in index[a[-k:]]:
            if a != b:
                olen = overlap(a, b, k)
                if olen > best_olen:
                    reada, readb = a, b
                    best_olen = olen
    return reada, readb, best_olen

"""
This function implements the greedy shortest common superstring algorithm using an accelerated version of pick_maximal_overlap(reads, k).
"""
def greedy_scs_index(reads, k):
    read_a, read_b, olen = pick_maximal_overlap_index(reads, k)
    while olen > 0:
        print(len(reads))
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = pick_maximal_overlap_index(reads, k)
    return ''.join(reads)