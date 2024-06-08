from Bio import SeqIO
import os


path = "/bulk/sycuro_bulk/daniel/diabetes/blast_picrust2"
os.chdir(path)
for file in os.listdir(path):
    if file.startswith('seq'):
        with open(file,'r') as open_file:
            out_file = open(path+"/just_one_seq/"+file, "w")
            for record in SeqIO.parse(open_file, "fasta"):
                out_file.write(f">{record.id}\n{record.seq}\n")
                break
