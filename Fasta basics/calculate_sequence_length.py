from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description = "This program will tell you the length of each fasta sequence found in a fasta file")
parser.add_argument("--file", type=str, help="The name of your fasta file")
args = parser.parse_args()
file = args.file
print(f"{file}")
for record in SeqIO.parse(file,"fasta"):
        print(f"{record.id}\t{len(record.seq)}")
