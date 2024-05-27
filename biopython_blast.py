from Bio.Blast import NCBIWWW, NCBIXML

# Define the unknown DNA sequence
unknown_dna = """TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTAC
AATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCAC
CTACGGTAGAG"""

# Perform a BLAST search against the nt (nucleotide) database
result_handle = NCBIWWW.qblast("blastn", "nt", unknown_dna)

# Parse the BLAST results
blast_records = NCBIXML.read(result_handle)

# Print the top BLAST hit
for alignment in blast_records.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < 0.01:  # Consider hits with e-value < 0.01
            print(f"****Alignment****")
            print(f"sequence: {alignment.title}")
            print(f"length: {alignment.length}")
            print(f"e-value: {hsp.expect}")
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
            break
    break