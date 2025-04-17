#Daniel Castaneda Mogollon, PhD
#April 17th, 2025
#This script renames .csv files by removing the word 'fasta' from it

import os
import pandas as pd

path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/t1d_db_fixed_discussed/fixed_assembly_benchling"
os.chdir(path)

for file in os.listdir():
    if '.fasta' in file:
        new_file = file.replace('.fasta', '')
        os.rename(file, new_file)
