#Daniel Castaneda Mogollon, PhD
#04/20/2025
#This script renames the .csv files and the .fasta files from the Sanger sequencing full-length and the
#blast results in .csv output. This matches the original sample names from UofT

import os
import pandas as pd

path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/t1d_db_fixed_discussed/fixed_assembly_benchling"
path_sanger = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/t1d_db_fixed_discussed/fixed_assembly_benchling/full_length"
os.chdir(path)
master_file = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/t1d_db_fixed_discussed/FINAL_report_WGS_all_csv_April.csv"
df = pd.read_csv(master_file)
subset_df = df[df['Selected for Downstream']=='Yes']
subset_df = subset_df[['Sample ID','Selected for Downstream']]
sample_id_list = subset_df['Sample ID'].apply(lambda x : x.split('_')[-1]).tolist() #Splits the sample IDs and takes only the number
my_dict = dict(zip(subset_df['Sample ID'], sample_id_list)) #Original sample ID + the number within sample ID, i.e. 007
for files in os.listdir(path):
    if files.endswith(".csv"):
        name = files.split('.')[0]
        name = name.split('seq')[1]
        for item in my_dict.items():
            if item[1] == name:
                new_name = item[0]+'.csv'
                #print('Old file '+str(files)+' and new file: '+str(new_name))
                os.renames(files, new_name)
os.chdir(path_sanger)
for files in os.listdir(path_sanger):
    if files!='.DS_Store':
        name = files.split('.')[0]
        name = name.split('seq')[1]
        for item in my_dict.items():
            if item[1] == name:
                new_name2 = item[0]+'.fasta'
                os.renames(files, new_name2)
                #print('Old fasta file name: '+files+' and new name: '+new_name2)
