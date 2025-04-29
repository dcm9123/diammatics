import pandas as pd
import os

path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/t1d_db_fixed_discussed/good_files_to_use/"
os.chdir(path)
master_file = "FINAL_report_WGS_r220_May_m.csv"
full_species_dict = {} #Making a dictionary to store values
master_df = pd.read_csv(master_file)
subset_master_df = master_df[["Sample ID","Classification","Selected for Downstream"]] #Subsetting only 3 columns
filtered_df = subset_master_df[subset_master_df["Selected for Downstream"]=="Yes"].copy() #Makes a new independent df in memory, so pandas knows what to change
filtered_df["Classification"] = filtered_df["Classification"].str.replace("d__","",regex=False)
filtered_df["Classification"] = filtered_df["Classification"].str.replace("p__","",regex=False)
filtered_df["Classification"] = filtered_df["Classification"].str.replace("c__","",regex=False)
filtered_df["Classification"] = filtered_df["Classification"].str.replace("o__","",regex=False)
filtered_df["Classification"] = filtered_df["Classification"].str.replace("f__","",regex=False)
filtered_df["Classification"] = filtered_df["Classification"].str.replace("g__","",regex=False)
filtered_df["Classification"] = filtered_df["Classification"].str.replace("s__","",regex=False)
full_species_dict = dict(zip(filtered_df["Sample ID"],filtered_df["Classification"]))
in_fasta = "t1d_db_m.fasta"
taxo_db = open("t1d_db_fullTaxo.fasta","w")
species_db = open("t1d_db_species.fasta","w")
with (open(in_fasta,"r")) as f:
    for line in f:
        line = line.replace("\n","")
        line = line.strip()
        if line.startswith(">"):
            line = line.replace(">","")
            line2 = line.replace("_v1v9","")
            line2 = line2.replace("_v3v4","")
            line2 = line2.replace("_069_2","_069")
            if line2 in full_species_dict:
                new_header = ">"+full_species_dict[line2].replace(" ","_")+"("+line+")"
                new_header2 = ">"+line+" "+full_species_dict[line2].split(";")[-1]
                print(new_header2)
            if line2 not in full_species_dict:
                print(f"Not found in full_species_dict: {line2}")
            #taxo_db.write(new_header+"\n")
            species_db.write(new_header2+"\n")
        else:
            #taxo_db.write(line+"\n")
            species_db.write(line+"\n")

