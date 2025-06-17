#Daniel Castaneda Mogollon, PhD
#June 17th, 2025
#This script simply merges two .tsv files by using the same column as overlap

import pandas as pd
import os

path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/"
os.chdir(path)
file1 = path+"t1d_db_fixed_discussed/FINAL_report_042525_r220.csv"
file2 = path+"Depth_of_coverage_genomes_Jun2025.csv"
df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

df_merged = pd.merge(left=df1, right=df2, how="outer", left_on="Sample ID", right_on="Sample_ID")
df_merged.to_csv("FINAL_report_061725_r220.csv")
