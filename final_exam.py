from Bio.Blast import NCBIWWW, NCBIXML

CODON_LEN=3

# Create dictionary of DNA sequences
def dna_sequences(filename):

    # Open fast examples
    seqs={}
    try:
        f=open(filename, 'r')
        name=''
        for line in f:
            line=line.rstrip()
            if line[0] == '>':
                words=line.split()
                name=words[0][1:] # remove the first character
                seqs[name]=''
            else:
                seqs[name]=seqs[name]+line
        return seqs
        f.close()
    except IOError:
        print('File doesnt exist')

# dna_seqs = dna_sequences('./dna.example_2.fasta')
dna_seqs = dna_sequences('./dna2.fasta')


"""
(1) How many records are in the file? A record in a FASTA file is defined as a single-line header, 
followed by lines of sequence data. The header line is distinguished from the sequence data by a greater-than (">")
symbol in the first column. The word following the ">" symbol is the identifier of the sequence,
and the rest of the line is an optional description of the entry. There should be no space between the ">"
and the first letter of the identifier. 
"""
number_of_records = len(dna_seqs)
# print(number_of_records)


"""
(2) What are the lengths of the sequences in the file? What is the longest sequence and what is the
shortest sequence? Is there more than onelongest or shortest sequence? What are their identifiers? 
"""

# Initialize variables to track the keys with the longest and shortest values
max_key = None
min_key = None
max_length = -1
min_length = float('inf')

# Loop through the dictionary and find the keys with the longest and shortest values
for key, value in dna_seqs.items():
    value_length = len(value)
    
    if value_length > max_length:
        max_length = value_length
        max_key = key
    
    if value_length < min_length:
        min_length = value_length
        min_key = key

# print(max_key, max_length)
# print(min_key, min_length)


"""
(3) In molecular biology, a reading frame is a way of dividing the DNA sequence of nucleotides
into a set of consecutive, non-overlapping triplets (or codons). Depending on where we start,
there are six possible reading frames: three in the forward (5' to 3') direction and
three in the reverse (3' to 5').

Given an input reading frame on the forward strand (1, 2, or 3) your program should be able to
identify all ORFs present in each sequence of the FASTA file, and answer the following questions: 
 - What is the length of the longest ORF in the file?
 - What is the identifier of the sequence containing the longest ORF?
 - For a given sequence identifier, what is the longest ORF contained in the sequence represented by that identifier?
 - What is the starting position of the longest ORF in the sequence that contains it?
    (i.e., The position should indicate the character number in the sequence)
"""

def find_orfs(sequence, frame):
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []
    
    # Adjust the sequence for the specified frame
    sequence = sequence[frame-1:]
    
    for i in range(0, len(sequence) - 2, CODON_LEN):
        codon = sequence[i:i+CODON_LEN]
        if codon == start_codon:
            for j in range(i+CODON_LEN, len(sequence) - 2, CODON_LEN):
                stop_codon = sequence[j:j+CODON_LEN]
                if stop_codon in stop_codons:
                    orf_position = i+frame
                    orfs.append((orf_position, sequence[i:j+CODON_LEN]))
                    break
    return orfs

orfs_dict = {}
FRAME = 1

# Loop through the dictionary and find the keys with the longest and shortest values
for key, value in dna_seqs.items():
    # if key == "gi|142022655|gb|EQ086233.1|16":
    #     orfs_dict[key] = find_orfs(value, FRAME)
    orfs_dict[key] = find_orfs(value, FRAME)

def find_longest_shortest_orf(dictionary):
    longest_position = None
    shortest_position = None
    longest_key = None
    shortest_key = None
    max_length = 0
    min_length = float('inf')

    for key, orfs in dictionary.items():
        for start_pos, orf_sequence in orfs:
            orf_length = len(orf_sequence)
            if orf_length > max_length:
                max_length = orf_length
                longest_key = key
                longest_position = start_pos
            if orf_length < min_length:
                min_length = orf_length
                shortest_key = key
                shortest_position = start_pos

    return longest_key, shortest_key, longest_position, shortest_position

longest_key, shortest_key, longest_position, shortest_position = find_longest_shortest_orf(orfs_dict)

# print(longest_position)
# print(shortest_position)
# print(longest_key)
# print(f"The key with the longest ORF is '{longest_key}' with a length of {max([len(orf[1]) for orf in orfs_dict[longest_key]])}.")
# print(f"The key with the shortest ORF is '{shortest_key}' with a length of {min([len(orf[1]) for orf in orfs_dict[shortest_key]])}.")

"""
(4) A repeat is a substring of a DNA sequence that occurs in multiple copies (more than one) 
somewhere in the sequence. Although repeats can occur on both the forward and reverse strands 
of the DNA sequence, we will only consider repeats on the forward strand here. Also we will 
allow repeats to overlap themselves. For example, the sequence ACACA contains two copies of 
the sequence ACA - once at position 1 (index 0 in Python), and once at position 3. 
Given a length n, your program should be able to identify all repeats of length n in 
all sequences in the FASTA file. Your program should also determine how many times 
each repeat occurs in the file, and which is the most frequent repeat of a given length.
"""


def find_most_frequent_repeat(sequences, n):
    all_repeats = {}
    for key, value in sequences.items():
        repeats = {}
        for i in range(len(value) - n + 1):
            repeat = value[i:i+n]
            if repeat in repeats:
                repeats[repeat] += 1
            else:
                repeats[repeat] = 1
        all_repeats[key] = repeats
    
    # Combine repeats from all sequences
    combined_repeats = {}
    for repeats in all_repeats.values():
        for repeat, count in repeats.items():
            if repeat in combined_repeats:
                combined_repeats[repeat] += count
            else:
                combined_repeats[repeat] = count
    
    # Find the most frequent repeat
    most_frequent_repeat = max(combined_repeats, key=combined_repeats.get)
    frequency = combined_repeats[most_frequent_repeat]
    
    return most_frequent_repeat, frequency

def count_max_repeats(max_count):
    all_repeats = {}
    
    for key, value in sequences.items():
        repeats = {}
        for i in range(len(value) - n + 1):
            repeat = value[i:i+n]
            if repeat in repeats:
                repeats[repeat] += 1
            else:
                repeats[repeat] = 1
        all_repeats[key] = repeats
    
    max_repeats = set()
    for repeats in all_repeats.values():
        for repeat, count in repeats.items():
            if count == max_count:
                max_repeats.add(repeat)
    return len(max_repeats)

# Example usage:
# Assuming 'sequences' is a dictionary where keys are sequence identifiers
# and values are DNA sequences.

# Find the most frequent repeat of length 6 in all sequences
most_frequent_repeat, frequency = find_most_frequent_repeat(dna_seqs, 7)

print(f"The most frequent repeat of length 7 is '{most_frequent_repeat}' "
      f"which occurs {frequency} times in all sequences.")


def find_repeats(sequences, n):
    all_repeats = {}
    for key, value in sequences.items():
        repeats = {}
        for i in range(len(value) - n + 1):
            repeat = value[i:i+n]
            if repeat in repeats:
                repeats[repeat] += 1
            else:
                repeats[repeat] = 1
        all_repeats[key] = repeats
    return all_repeats

def count_max_repeats(all_repeats, max_count):
    max_repeats = set()
    for repeats in all_repeats.values():
        for repeat, count in repeats.items():
            if count == max_count:
                max_repeats.add(repeat)
    return len(max_repeats)

# Example usage:
# Assuming 'sequences' is a dictionary where keys are sequence identifiers
# and values are DNA sequences.

# Find all repeats of length 12 in the input file
all_repeats = find_repeats(dna_seqs, 12)

# Find the maximum count of the most frequent repeat of length 12
max_count = max(max(repeats.values(), default=0) for repeats in all_repeats.values())

# Count the number of different 12-base sequences that occur 'max_count' times
num_different_sequences = count_max_repeats(all_repeats, max_count)

print(f"The number of different 12-base sequences that occur {max_count} times is: {num_different_sequences}")
