
from collections import defaultdict

def overlap(a, b, min_length=3):
    """Return length of longest suffix of 'a' matching
    a prefix of 'b' that is at least 'min_length' characters long.
    If no such overlap exists, return 0."""
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def read_fastq_sequences(filename):
    """Read sequences from a FASTQ file (assumes every 4th line is a sequence)."""
    sequences = []
    with open(r"C:\Users\liamm\Desktop\ERR266411_1.for_asm.fastq") as f:
        while True:
            f.readline()  # skip header
            seq = f.readline().strip()
            if not seq:
                break
            sequences.append(seq)
            f.readline()  # skip +
            f.readline()  # skip quality
    return sequences

def build_kmer_index(reads, k):
    """Build an index mapping each k-mer to the set of reads containing it."""
    index = defaultdict(set)
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            index[kmer].add(read)
    return index

def count_overlaps(reads, k=30):
    """Count the number of read pairs with an exact suffix-prefix overlap of at least k bases."""
    kmer_index = build_kmer_index(reads, k)
    overlaps = set()

    for read in reads:
        suffix = read[-k:]
        for candidate in kmer_index[suffix]:
            if read != candidate:
                olen = overlap(read, candidate, min_length=k)
                if olen > 0:
                    overlaps.add((read, candidate))

    return len(overlaps)

# Usage
filename = "ERR266411_1.for_asm.fastq"
reads = read_fastq_sequences(filename)
overlap_count = count_overlaps(reads, k=30)
print(f"Total number of overlaps (edges in the graph): {overlap_count}")

def count_overlaps_and_sources(reads, k=30):
    """Count total overlaps and number of reads with at least one outgoing edge."""
    kmer_index = build_kmer_index(reads, k)
    overlaps = set()
    source_reads = set()

    for read in reads:
        suffix = read[-k:]
        for candidate in kmer_index[suffix]:
            if read != candidate:
                olen = overlap(read, candidate, min_length=k)
                if olen > 0:
                    overlaps.add((read, candidate))
                    source_reads.add(read)  # This read has an outgoing edge

    return len(overlaps), len(source_reads)

filename = "ERR266411_1.for_asm.fastq"
reads = read_fastq_sequences(filename)
total_edges, total_sources = count_overlaps_and_sources(reads, k=30)
print(f"Total overlaps (edges): {total_edges}")
print(f"Reads with at least one outgoing edge: {total_sources}")
