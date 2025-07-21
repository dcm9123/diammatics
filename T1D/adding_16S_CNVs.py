# Daniel Castaneda Mogollon, PhD
# 21:23
# adding_16S_CNVs.py.py
# Purpose: This script takes the original 116 samples from the MASTER file and creates a new 16S.txt file for
# PICRUSt2. Because my 16S database has ~209 samples (Sanger, in silico, barrnap, etc), there will be 'repetitions' of 16S CNV
# that must be added for PICRUSt2 to work

import os
import pandas as pd
from Bio import SeqIO
import re

#Getting the MASTER file with the final 16S count:
path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/"
os.chdir(path)
input_file = pd.read_csv("FINAL_report_062325_r220_CNVs.csv")
df = input_file[["ID","Headers","Final 16S number"]]
#df = df.set_index("ID")


#Getting the 16S db file with the 209 IDs
reference_16_list = []
ID_file = os.path.join(path+"picrust2_june232025/reference_16S.fasta")
with open(ID_file,'r') as handle:
    for record in SeqIO.parse(handle, "fasta"):
        reference_16_list.append(record.id)

#Matching the ID fasta with the MASTER file 16S numbers
ids_dict = {}
for name in reference_16_list:
    name_m = re.split('_',name)[0:4]
    name_m = "_".join(name_m)
    ids_dict[name] = name_m

i = 0
f_out = open("picrust2_june232025/16S.txt","w")
f_out.write("assembly\t16S_rRNA_Count\n")
id_list = df["ID"].tolist()
for key,value in ids_dict.items():
    if value in id_list:
        index_val = df[df["ID"] == value].index[0] #Retrieves the corresponding index where the value matches the excel file
        val_16S = df.loc[index_val,"Final 16S number"] #Gets the corresponding 16s NUMBER
        f_out.write(f"{key}\t{val_16S}\n")
