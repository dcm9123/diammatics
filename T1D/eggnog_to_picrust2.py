# Daniel Castaneda Mogollon, PhD
# 14:49
# eggnog_to_picrust.py.py
# Purpose: This script was made to format the output from emapper to something picrust2 can understand.
#It is inteded to format EC and KOs.
import re
from re import split
import pandas as pd
import os


#This allows me to print the entire data frame if needed
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 0)  # Let pandas auto-detect width

#Calling the first function to format the original eggnog file to something pandas-readable
def formatting_original_file(path):
    os.chdir(path)
    for file in os.listdir():
        if file.endswith(".annotations"):
            file_name = file.split(".")[0]
            f_out = "tmp.txt"
            top_to_remove = 4 #Removes the first 4 lines of the eggnog emapper file
            last_to_remove = 3 #Removes the last 3 lines of the eggnog emapper file
            with open(file,"r") as infile:
                lines = infile.readlines()
            trimmed_lines = lines[top_to_remove:]
            trimmed_lines = trimmed_lines[:-last_to_remove]
            with open(f_out,"w") as outfile:
                outfile.writelines(trimmed_lines)           #Writes the new file without the useless lines
                # Reading the output file as a df now
                df = pd.read_table(f_out)
                df.insert(0, "Sample ID", file_name)  # Inserts the new column in the very first position
                df.to_csv(file_name + "_eggnog.tsv", sep='\t', index=False)  # Output as tsv file
            os.remove(f_out)
        else:
            print(file+" is not an emapper.annotations file, skipping...")

def formatting_picrust2_annotations(path):
    os.chdir(path)
    file = "S_S5_Blss_094_eggnog.tsv"
    prefix = file.split("_eggnog")[0] #Removes the _eggnog.tsv as part of the variable. This will be used to name the genomes
    df = pd.read_csv(file,sep='\t') #Reading the file
    df["EC"] = df["EC"].replace("-","NA")
    df["EC"] = df["EC"].astype(str).str.replace(r'(\d+(?:\.\d+)+)', r'EC:\1', regex=True)
    df["KEGG_ko"] = df["KEGG_ko"].replace("-","NA") #Replaces entire cell containing '-'
    df["KEGG_ko"] = df["KEGG_ko"].str.replace("ko:","") #Replaces instances of "ko:" in a cell, that's why we used str.replace
    df_subset = df.iloc[:,[0,11,12]].copy() #This calls only the ID of the sample, the ECs, and the KOs, and makes sure I work on a copy and do not modify the original df
    print(df_subset)
    for column in df_subset.columns[1:3]:
        df_subset[column] = df_subset[column].astype(str) #Makes sure that the EC numbers are treated as a string
        df_subset[column] = df_subset[column].str.split(",") #Separates the multiple EC numbers of each cell into different rows
        df_exploded_ec_ko = df_subset.explode(column) #Separates into new rows the split values
        ec_ko_counts = df_exploded_ec_ko[column].value_counts() #Counts the instances of each EC value
        ec_ko_new = pd.DataFrame([ec_ko_counts.to_dict()]) #Turns these into a dictionary, where the key is the EC number, and the value is the count. It is transformed from series to df
        ec_ko_new.drop(columns="NA", inplace=True) #Deletes the 'NA' value column
        ec_ko_new.insert(0,"Sample ID",prefix)
        print(ec_ko_new)

def main():
    formatting_original_file("/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2_june232025/genome_species/eggnogs/annotations")
    formatting_picrust2_annotations("/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2_june232025/genome_species/eggnogs/annotations")
if __name__ == "__main__":
    main()
