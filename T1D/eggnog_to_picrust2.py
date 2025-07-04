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
        file_name = file.split(".")[0]
        f_out = "tmp.txt"
        top_to_remove = 4
        last_to_remove = 3
        with open(file,"r") as infile:
            lines = infile.readlines()
        trimmed_lines = lines[top_to_remove:]
        trimmed_lines = trimmed_lines[:-last_to_remove]
        with open(f_out,"w") as outfile:
            outfile.writelines(trimmed_lines)

    #Reading the output file as a df now
        df = pd.read_table(f_out)
        df.insert(0,"Sample ID",file_name) #Inserts the new column in the very first position
        df.to_csv(file_name+"_eggnog.tsv",sep='\t',index=False) #Output as tsv file
    os.remove(f_out)
path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2_june232025/genome_species/eggnogs/annotations"
def main():
    formatting_original_file(path)
if __name__ == "__main__":
    main()
