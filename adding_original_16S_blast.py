import pandas as pd
import os
from Bio import SeqIO

path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/t1d_db_fixed_discussed"
os.chdir(path)
dict_ids = {}
for file in os.listdir():
    if (file.endswith(".csv") and ("FINAL" or "Final") not in file):
        print(file)
        df = pd.read_csv(file)
        for seq_id in df["Sequence ID"]:
            dict_ids[seq_id] = 0
        name = file.split("blast")[0]
        for file2 in os.listdir(path+"/16S_from_genomes"):
            if name in file2:
                in_fasta = SeqIO.parse(path+"/16S_from_genomes/"+file2, "fasta")
                for record in in_fasta:
                    if record.id in dict_ids:
                        dict_ids[record.id] = str(record.seq)
        df['Original Sequence'] = df["Sequence ID"].map(dict_ids)
        df.to_csv(name+"_blast_filtered_with_original_seq.csv")


