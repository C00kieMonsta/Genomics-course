from Bio.Seq import Seq

# Define the DNA sequence
dna_sequence = """TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTAC
AATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCAC
CTACGGTAGAG"""

# Remove any newlines or spaces in the sequence
dna_sequence = dna_sequence.replace('\n', '').replace(' ', '')

# Create a Seq object
my_seq = Seq(dna_sequence)

# Translate the DNA sequence into a protein sequence
protein_sequence = my_seq.translate()

# Print the protein sequence
print(protein_sequence)
