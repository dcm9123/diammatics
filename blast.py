from Bio.Blast import NCBIWWW
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="This program was created with the purpose of running blast on a multifasta file")
parser.add_argument("--database", type=str, help="The name of the database that is taken from the NCBI website")
parser.add_argument("--file", type=str, help="This is the name and path of your multifasta file")

args = parser.parse_args()
database = args.database
fasta_file = args.file

def blast_fasta(fasta_file, database):
    sequences = list(SeqIO.parse(fasta_file,"fasta"))
    for seq in sequences:
        print(seq.seq)
        result_handle = NCBIWWW.qblast("blastn", database, seq.seq)
        output_file = seq.id + "_blast.xml"
        with open(output_file,"w") as out_handle:
            out_handle.write(result_handle.read())
        result_handle.close()

        print(f"The results for {seq.id} saved to {output_file}")

blast_fasta(fasta_file,database)
