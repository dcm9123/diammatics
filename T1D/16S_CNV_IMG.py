# Daniel Castaneda Mogollon, PhD
# 11:04
# 16S_CNV_IMG.py
# Purpose: This script was made with the purpose of extracting basic 16S rRNA information from the rrnDB database.
# This data will be used for the T1D PICRUSt2 functional analyses.

import pandas as pd
import os

#This allows me to print the entire df without any issues or constraints from python
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 0)  # Let pandas auto-detect width

path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes"
os.chdir(path)
#In here, I am reading the rrnDB file with the 16S information of all the species in their latest version
df = pd.read_csv("rrnDB-5.10.csv", dtype={"Data source" : "string", "NCBI tax id":"Int64", "BioSample":"string",
                                          "BioProject":"string","Data Source":"string","NCBI scientific name":"string",
                                          "RDP taxa":"string","RDP taxonomy":"string","basecount":"Int64","16S gene count":"Int64",
                                          "23S gene count":"Int64","tRNA gene count":"Int64","Evidence":"string","Note":"string",
                                          "References":"string"})
#Reading my own master file
df_master = pd.read_csv(path+"/t1d_db_fixed_discussed/good_files_to_use/FINAL_report_WGS_r220_June.csv")
df_depth = pd.read_csv(path+"/t1d_db_fixed_discussed/good_files_to_use/depth_coverage.csv")

#This are all the species I have in my master file across all 4 consortia
species_list = ["Agathobacter faecis","Akkermansia massiliensis","Akkermansia muciniphila","Alistipes finegoldii",
                "Alistipes onderdonkii","Alitiscatomonas sp900066535","Lacrimispora amygdalina","Bacillus subtilis",
                "Bacteroides stercoris","Bacteroides thetaiotaomicron","Bacteroides uniformis","Bacteroides xylanisolvens",
                "Bifidobacterium adolescentis","Bifidobacterium dentium","Bifidobacterium longum","Bifidobacterium pseudocatenulatum",
                "CHH4-2 sp018378255","Collinsella sp902362275","Eggerthella lenta","Eisenbergiella porci","Eisenbergiella tayi",
                "Enterocloster alcoholdehydrogenati","Enterocloster bolteae","Enterocloster citroniae","Enterocloster clostridioformis",
                "Enterocloster sp005845215","Enterococcus faecalis","Enterococcus_A avium","Enterococcus_B durans","Escherichia coli",
                "Faecalibacterium longum","Flavonifractor plautii","Holdemania filiformis","Hungatella effluvi","Intestinomonas_A celer",
                "Otoolea symbiosa","Parabacteroides distasonis","Parabacteroides merdae","Phocaeicola dorei","Phocaeicola vulgatus",
                "Pseudoflavonifractor_A sp022772585","Roseburia_C amylophila","Ruminococcus_B gnavus","Sarcina perfringens",
                "Staphylococcus capitis","Staphylococcus epidermidis","Staphylococcus hominis","Sterptococcus sp001556435",
                "Sutterella wadsworthensis"]

#This are alternative names for some of the same species by using old nomenclature or NCBI nomenclature instead of GTDB
species_list_alt = ["Lacrimospora amygdalina","Clostridium sp. Chh4-2","Collinsella aerofaciens","Enterococcus avium","Enterococcus durans",
                    "Lawsonibacter celer","Clostridium symbiosum","Pseudoflavonifractor gallinarum","Roseburia amylophila",
                    "Ruminococcus gnavus","Clostridium perfringens","Roseburia faecis"]

#Filtering names by the actual genomes we used
df_master_filtered = df_master[df_master['Selected for Downstream']=="Yes"]
species_list_master = df_master_filtered['GTDBtk Species Classification'].unique().tolist() #Counting the number of unique species (removes repeated names)
species_overlap = []
species_alt_overlap = []

i=0
#print("There are "+str(len(species_list_master))+" unique species in my original list, but the ones missing are: \n")
for item in species_list_master:
    if item not in df['NCBI scientific name'].tolist():
        i=i+1
    else:
        species_overlap.append(item) #If there's a match between the names of my master file and rrnDB, then they are retained and appended into a new list

i=0
for item in species_list_alt:
    if item not in df['NCBI scientific name'].tolist(): #Same but for the alternative list
        i=i+1
        #print(item)
    else:
        species_alt_overlap.append(item)
#print("There are "+str(i)+"/"+str(len(species_list_alt))+" species with alternative names missing from rrnDB in the NCBI list:")

list_of_dicts=[]
for list in [species_overlap,species_alt_overlap]: #Iterating over the two lists I have with overlap in names
    overlap_dict={}
    for item in list: #This section just consideres the NCBI scientific name, and 16S info from the rrnDB. I don't need to groupby species as I am already filtering
        df_filtered = df[df['NCBI scientific name']==item]
        data_16S = df_filtered['16S gene count']
        median_16S = data_16S.median() #Calculation of the same species median value and the other values.
        min_16S = data_16S.min()
        max_16S = data_16S.max()
        count_16S = data_16S.count()
        overlap_dict[item] = {'Count':count_16S,'Median 16S':median_16S,'Minimum 16S':min_16S,'Maximum 16S':max_16S} #Appending the dictionary with keys and values
    df_16S = pd.DataFrame.from_dict(overlap_dict, orient='index').reset_index() #Converting the dictionary into a data frame
    df_16S.rename(columns={'index':'Species rrnDB','Count':'16S count'}, inplace=True) #Renaming columns
    list_of_dicts.append(df_16S) #Appending df1 and df2 for each list

df_rrndb = pd.concat(list_of_dicts).reset_index(drop=True) #This joins the two data frames (or stacks them one on top of the other) and resets the index
df_rrndb.to_csv("rrndb_16S.csv",index=False) #Writing final df of info needed

#Mergind the df I just made with the master file only if the same species overlap between the two files.
df_joined = pd.merge(df_master_filtered,list_of_dicts[0],how="outer",left_on="GTDBtk Species Classification",right_on="Species rrnDB").reset_index(drop=True)
df_joined = pd.merge(df_joined, list_of_dicts[1],how="outer",left_on="Other potential names", right_on="Species rrnDB").reset_index(drop=True)
df_joined = pd.merge(df_joined, df_depth,how="outer",left_on="Sample ID",right_on="ID")

#I had to manually delete the extra column Species rrnDB_x and 'y' into one, as I have two different species names. I also added 4 columns for IMG data
df_joined.set_index("Sample ID",inplace=True)
df_joined.to_csv("test16S.csv",index=False)

#df1_high_q = df1[(df1['CheckM2 Contamination']<5.0) & (df1['CheckM2 Completeness']>95.0) & (df1['High Quality']=="Yes")]
#df1_subset = df1_high_q[['Genome Name / Sample Name', 'IMG Genome ID', 'NCBI Species', 'GTDB Species', 'Assembly Method',
#                  'High Quality', '16S rRNA Count  * assembled', 'CheckM2 Completeness', 'CheckM2 Contamination']]
##df1_subset = df1.groupby('NCBI Species')['16S rRNA Count  * assembled'].agg(['count','mean','median','min','max']).reset_index()
#df1_subset_final = df1.merge(df1_subset, on = 'NCBI Species')
#df1_final = df1_subset_final[['IMG Genome ID','Genome Name / Sample Name','NCBI Species','GTDB Species','High Quality','CheckM2 Completeness','CheckM2 Contamination','count','mean','median','min','max']]
