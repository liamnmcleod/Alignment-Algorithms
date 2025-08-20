"""Script to find a reverse complement of DNA, then apply a naive matching algorithm + reverse complement and up
to 2 mismatch positions"""

#returns reverse complement of DNA
def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


# naive exact matching algorithm
# Returns a list of occurrences (offset)
def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

# naive exact matching algorithm that is strand-aware
# Instead of looking only for occurrences of P in T, additionally look for occurrences of the reverse complement of P in T.
def naive_with_rc(p, t):
	r = reverseComplement(p)
	if r == p:
		return naive(p,t)
	else:
		return naive(p,t) + naive(r,t)

# naive exact matching algorithm
# Returns a list of occurrences (offset)
def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        mismatch_count = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mismatch_count+=1

                if mismatch_count > 2:
                	match = False
                	break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

# parses a DNA reference genome from a file in the FASTA format.
def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

# parses the read and quality strings from a FASTQ file containing sequencing reads
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

p = 'CCC'
ten_as = 'AAAAAAAAAA'
t = ten_as + 'CCC' + ten_as + 'GGG' + ten_as
occurrences = naive_with_rc(p, t)
print(p)
print(t)
print(occurrences)
# Assert: [10, 23]


p = 'CGCG'
t = ten_as + 'CGCG' + ten_as + 'CGCG' + ten_as
occurrences = naive_with_rc(p, t)
print(p)
print(t)
print(occurrences)
# Assert: [10, 24]


# wget http://d396qusza40orc.cloudfront.net/ads1/data/phix.fa
phix_genome = readGenome('phix.fa')
occurrences = naive_with_rc('ATTA', phix_genome)
print('offset of leftmost occurrence: %d' % min(occurrences))
# offset of leftmost occurrence: 62

print('# occurrences: %d' % len(occurrences))
# occurrences: 60

# Question 1
# How many times does AGGT or its reverse complement (ACCT) occur in the lambda virus genome?
# E.g. if AGGT occurs 10 times and ACCT occurs 12 times, you should report 22.
lambda_virus_genome = readGenome('lambda_virus.fa')
occurrences = naive_with_rc('AGGT', lambda_virus_genome)
print('AGGT in lambda_virus:')
print('# occurrences: %d' % len(occurrences))

# Question 2
# How many times does TTAA or its reverse complement occur in the lambda virus genome?
# Hint: TTAA and its reverse complement are equal, so remember not to double count.
occurrences = naive_with_rc('TTAA', lambda_virus_genome)
print('TTAA in lambda_virus:')
print('# occurrences: %d' % len(occurrences))


# Question 3
# What is the offset of the leftmost occurrence of ACTAAGT or its reverse complement in the Lambda virus genome?
# E.g. if the leftmost occurrence of ACTAAGT is at offset 40 (0-based) and the leftmost occurrence of its reverse
# complement ACTTAGT is at offset 29, then report 29.
occurrences = naive_with_rc('ACTAAGT', lambda_virus_genome)
print('offset of leftmost occurrence of ACTAAGT: %d' % min(occurrences))

# Question 4
# What is the offset of the leftmost occurrence of AGTCGA or its reverse complement in the Lambda virus genome?
occurrences = naive_with_rc('AGTCGA', lambda_virus_genome)
print('offset of leftmost occurrence of AGTCGA: %d' % min(occurrences))

# Question 5
# As we will discuss, sometimes we would like to find approximate matches for P in T. T
# hat is, we want to find occurrences with one or more differences.
# For Questions 5 and 6, make a new version of the naive function called naive_2mm that allows up to
# 2 mismatches per occurrence. Unlike for the previous questions, do not consider the reverse complement here.
# We're looking for approximate matches for P itself, not its reverse complement. ￼
# For example, ACTTTA occurs twice in ACTTACTTGATAAAGT, once at offset 0 with 2 mismatches, and once at offset 4 with 1 mismatch.
# So naive_2mm('ACTTTA', 'ACTTACTTGATAAAGT') should return the list [0, 4].

# How many times does TTCAAGCC occur in the Lambda virus genome when allowing up to 2 mismatches?
occurrences = naive_2mm('TTCAAGCC', lambda_virus_genome)
print('TTCAAGCC (up to 2 mismatches) in lambda_virus:')
print('# occurrences: %d' % len(occurrences))

# Question 6
# What is the offset of the leftmost occurrence of AGGAGGTT in the Lambda virus genome when allowing up to 2 mismatches?
occurrences = naive_2mm('AGGAGGTT', lambda_virus_genome)
print('AGGAGGTT (up to 2 mismatches) in lambda_virus:')
print('offset of leftmost occurrence: %d' % min(occurrences))

# Test 1
p = 'CTGT'
ten_as = 'AAAAAAAAAA'
t = ten_as + 'CTGT' + ten_as + 'CTTT' + ten_as + 'CGGG' + ten_as
occurrences = naive_2mm(p, t)
print(occurrences)
# [10, 24, 38]


# Test 2
phix_genome = readGenome('phix.fa')
occurrences = naive_2mm('GATTACA', phix_genome)
print('offset of leftmost occurrence: %d' % min(occurrences))
# offset of leftmost occurrence: 10

print('# occurrences: %d' % len(occurrences))
# occurrences: 79






