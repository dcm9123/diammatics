# Daniel Castaneda Mogollon, PhD
# 21:54
# merging_functional_predictions.py.py
# Purpose: This script was generated to merge the results of the EC and KO functional prediction files from
# the 4 consortia (NS1, NS6, S2, S5) into one.

import os
import pandas as pd


#This allows me to print the entire df without any issues or constraints from python
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 0)  # Let pandas auto-detect width

path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2_june232025/"
os.chdir(path)
consortia_list = ["ns1","ns6","s2","s5"]
function_list = ["KO","EC"]
df_list_ko = []
df_list_ec = []
i = 0

for consortia in consortia_list:
    for function in function_list:
        df = pd.read_csv(f"{consortia}_output/{consortia}_{function}_metagenome_out/pred_metagenome_unstrat.tsv",sep='\t',header=0)
        df = df.set_index('function')
        if function == "KO":
            df_list_ko.append(df)
            i = i+1
        else:
            df_list_ec.append(df)
            i = i+1

df_ko_merged = pd.concat(df_list_ko, axis=1).fillna(0)
df_ec_merged = pd.concat(df_list_ec, axis=1).fillna(0)

df_ko_merged.to_csv("KO_merged_metagenome.tsv",index=True, sep="\t")
df_ec_merged.to_csv("EC_merged_metagenome.tsv",index=True, sep="\t")


for list_in in df_list_ko,df_list_ec:
    for dframe in list_in:
        print(f"The dimension of each df are: {dframe.shape}")

print(f"The dimensions of the merged KO df are {df_ko_merged.shape}")
print(f"The dimensions of the merged EC df are {df_ec_merged.shape}")



for item in df_ko_merged.index:
    found = any(item in df.index for df in df_list_ko)
    if not found:
        print(f"{item} not found in any KO input dataframe (unexpected)")

for dframe in df_list_ko:
    for value in dframe.index:
        if value not in df_ko_merged.index:
            print(f"{value} not found in the merged KO df (unexpected)")


for item in df_ec_merged.index:
    found = any(item in df.index for df in df_list_ec)
    if not found:
        print(f"{item} not found in any EC input dataframe (unexpected)")

for dframe in df_list_ec:
    for value in dframe.index:
        if value not in df_ec_merged.index:
            print(f"{value} not found in the merged EC df (unexpected)")