

def read_fasta():
    sequences = {}
    with open(r"C:\Users\liamm\Desktop\multifast.fasta", 'r') as fasta_file:
        for line in fasta_file: #reads line by line
            if line.startswith('>'): #check for lines starting with >, indicates seq_ID header
                seq_id = line.strip() #removes trailing whitespace
                seq_id = seq_id.replace('>', '') #removes > sign
            else:
                if seq_id not in sequences: #checks if the seq_id key is in the dictionary
                    sequences[seq_id] = ''
                sequences[seq_id] += line.strip()  #append to the dictionary: dict[key] = value
    return sequences #return the dictionary


sequences = read_fasta()
for seq_id in sequences:
    print(seq_id, (f":{sequences[seq_id]}"))

print("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
print(f"Number of sequences: {len(sequences)}")


#(Need to combine the lengths into the printed strings to make more pythonic)

#sequence_lengths = {seq_id: len(seq) for seq_id, seq in sequences.items()}
#for seq_id, length in sequence_lengths.items():
    #print(f"{seq_id} Sequence length: {length}")


def get_length(seq_id):
    return len(sequences[seq_id])
longest_seq_id = max(sequences)
shortest_seq_id = min(sequences)
print(f"Longest sequence length: {get_length(longest_seq_id)}")
print(f"Shortest sequence length: {get_length(shortest_seq_id)}")


