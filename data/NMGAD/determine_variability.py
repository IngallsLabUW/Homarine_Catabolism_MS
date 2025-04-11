from Bio import SeqIO, pairwise2
import numpy as np

# 1️⃣ Load sequences from a FASTA file
fasta_file = "/mnt/c/SynologyDrive/UW/UW_Anitra_Ingalls/Homarine_project/Data/NMGAD/NMGAD_SOX_genes/sox_genes/matchedsoxA_mgdC_gene_sequences.fasta"  # Change this to your filename
sequences = [str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]

print(f"Loaded {len(sequences)} sequences.")

# 2️⃣ Compute pairwise sequence identity
identities = []
num_seqs = len(sequences)

for i in range(num_seqs):
    for j in range(i + 1, num_seqs):
        alignments = pairwise2.align.globalxx(sequences[i], sequences[j], one_alignment_only=True)
        alignment = alignments[0]  # Take the best alignment
        matches = sum(a == b for a, b in zip(alignment.seqA, alignment.seqB))
        identity = (matches / len(alignment.seqA)) * 100
        identities.append(identity)

# Compute and print average identity
if identities:
    avg_identity = np.mean(identities)
    print(f"Average Pairwise Identity: {avg_identity:.2f}%")
else:
    print("Not enough sequences to compute pairwise identity.")
