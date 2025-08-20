#Naive exact algorithm for searching on both the forward and reverse complement
# r = read, t = genome template
from collections import Counter


filename = r"C:\Users\liamm\Desktop\chr1.GRCh38.excerpt.fasta"

def readGenome(filename):
    genome = ''
    with open(r"C:\Users\liamm\Desktop\chr1.GRCh38.excerpt.fasta", 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
print(f"Forward strand: {readGenome(filename)}")

#genome: str = readGenome(filename)


def reverseComplement(genome):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in genome:
        t = complement[base] + t
    return t
print(f"Reverse Strand: {reverseComplement(readGenome(filename))}")


def naive(r, t):
    occurrences = []
    for i in range(len(t) - len(r) + 1):  # loop over alignments
        match = True
        for j in range(len(r)):  # loop over characters
            if t[i+j] != r[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

matches = naive('ACTAAGT', )
print(f"Total matches: {len(matches)}")
print(naive('ACTAAGT', ''))


def reverse_naive(r, t):
    reverse_occurrences = []
    for i in range(len(t) - len(r) + 1):
        match = True
        for j in range(len(r)):
            if t[i+j] != r[j]:
                match = False
                break
        if match:  # This now runs once per alignment, not per base
            reverse_occurrences.append(i)
    return reverse_occurrences


reverse_matches = reverse_naive('ACTAAGT','')
print(f"Total reverse complement matches: {len(reverse_matches)}")
print(reverse_naive('ACTTAGT', ''))


def naive_2mm(p, t):
    mm_occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        mismatches = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:
                mismatches += 1
                if mismatches > 2:
                    break  # if too many mismatches stop
        if mismatches <= 2:
            mm_occurrences.append(i)  # accept this match
    return mm_occurrences

mm_matches = naive_2mm('TTCAAGCC','')
print(f"Total number of matches with up to 2mm: {len(mm_matches)}")
print(naive_2mm('TTCAAGCC',''))

