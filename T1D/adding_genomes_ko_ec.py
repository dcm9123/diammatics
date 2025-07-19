# Daniel Castaneda Mogollon, PhD
# 13:44
# adding_genomes_ko_ec.py.py
# Purpose: This script was made in order to add genomes from the same 16S to the KO and EC script. Because PICRUSt2
# needs an identical number of 16S and genomes annotated, I have to add 'repetitions' of these annotations based on the
# 16S sequence names (in my case I have multiple 16S from the same strains from Sanger or barrnap).

#Input: The user needs to provide a path where the reference 16S database is found (which should be a .fasta file called
# reference_16S.fasta")

import os

import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import Counter

def setting_up():
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 0)  # Let pandas auto-detect width
    return()

def listing_genomes(path):
    reference_16_list = []
    file = os.path.join(path,"reference_16S.fasta")
    with open(file,'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            reference_16_list.append(record.id)
    #print(reference_16_list)
    return(reference_16_list)


def getting_files(path):
    file_ko = pd.read_csv(os.path.join(path,'Final_KO_file.tsv'), sep='\t', header=0, encoding='utf-8')
    file_ec = pd.read_csv(os.path.join(path,'Final_EC_file.tsv'), sep='\t')
    file_ko.to_csv("test.tsv", sep='\t', index=False)
    file_ko_transposed = file_ko.transpose()
    #print(file_ko)
    genome_list = file_ko_transposed.iloc[0].tolist()
    return(genome_list,file_ko,file_ec)

def adding_names(reference_16S_list,genome_list, file_ko, file_ec):
    #reference_16S_list is the list of headers I have for my 16S database (209)
    #genome_list is the number of genome IDs found in the final KO fila (116)
    #file_ko has the entire 116 genomes annotated
    matches = []
    no_matches = []
    genome_dict = {}
    for item in reference_16S_list:
        if item in genome_list:
            genome_dict[item] = item
            matches.append(item)

        else:
            item_m = item.split("_")[:4]
            item_m = "_".join(item_m)
            genome_dict[item] = item_m
            no_matches.append(item)
    #In this part, I make sure that the genome dictionary has the same length as the reference_list
    #genome_dict keys are the headers of the 209 fasta sequences, and the values are the modified version that matches the KO and EC files
    counter1 = Counter(matches)
    counter2 = Counter(no_matches)
    for counter in counter1,counter2:
        for item,count in counter.items():
            if count > 1:
                print(f"{item}: {count}")
    print(f"There are a total of {len(genome_dict.values())} original header names in the dictionary values")
    print(f"There are a total of {len(file_ko)} IDs in the KO annotation file")

    #This will iterate over the EC file and the KO file
    for file_to_do,name in zip([file_ko,file_ec],["KO","EC"]): #Gives them appropriate names
        ko_indexed = file_to_do.set_index("assembly") #This makes the assemblies the actual index
        output_rows = []
        for header in genome_dict: #header is the actual key of the dictionary, the 209 sequences I have from my 16S db
            ko_id = genome_dict[header] #This is the actual name that will match the annotated genome file from eggnog
            if ko_id in ko_indexed.index:
                row = ko_indexed.loc[ko_id].copy() #Copying the row that belongs to the ko_id
            row.name = header #This names the assembly just like the original database
            output_rows.append(row) #appending into the list
        final_df = pd.DataFrame(output_rows) #Conveerting the list into a data frame
        final_df.index.name = "assembly" #Making sure the header says 'assembly'
        final_df = final_df.replace(r"^\s*$", 0, regex=True) #Replacing empty cells with 0
        final_df = final_df.fillna(0) #Replacing NAs with zeroes.
        final_df.to_csv(f"/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2_june232025/{name}_for_picrust2.tsv",sep='\t', index=True)

def main():
    setting_up()
    reference_16S_list = listing_genomes('/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2_june232025/')
    genome_list,file_ko,file_ec = getting_files('/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2_june232025/eggnogs_annotations/FINAL/')
    adding_names(reference_16S_list,genome_list, file_ko, file_ec)

if __name__ == '__main__':
    main()
