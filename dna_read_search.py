"""Simple script to locate substring matches in nucleic acid sequences."""

from Bio.Seq import Seq

filename = r"C:\Users\liamm\Desktop\lambda_virus.fa"

def read_genome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                genome += line.strip()
    return genome

# Read genome
dna = read_genome(filename)
dna_seq = Seq(dna)

# Define NGS read and its reverse complement
ngs_read = 'AGGT'
seq_ngs_read = Seq(ngs_read)
reverse_ngs_read = seq_ngs_read.reverse_complement()

# Search function that returns matches
def search_pattern_positions(pattern, text):
    m = len(pattern)
    n = len(text)
    positions = []
    for i in range(n - m + 1):
        if text[i:i+m] == pattern:
            positions.append(i)
    return positions

# Find positions
forward_positions = search_pattern_positions(str(seq_ngs_read), str(dna_seq))
reverse_positions = search_pattern_positions(str(reverse_ngs_read), str(dna_seq))

# Find the leftmost (smallest index) if any
all_positions = []
if forward_positions:
    all_positions.append(min(forward_positions))
if reverse_positions:
    all_positions.append(min(reverse_positions))

if all_positions:
    leftmost = min(all_positions)
    print(f"Leftmost occurrence at offset {leftmost} (0-based)")
else:
    print("No matches found for read or its reverse complement.")
