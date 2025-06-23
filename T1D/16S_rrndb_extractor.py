# Daniel Castaneda Mogollon, PhD
# 11:04
# 16S_CNV_IMG.py
# Purpose: This script was made with the purpose of extracting basic 16S rRNA information from the IMG database.
# This data will be used for the T1D PICRUSt2 functional analyses.

import pandas as pd
import os

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 0)  # Let pandas auto-detect width

path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes"
os.chdir(path)
df = pd.read_csv("rrnDB-5.10.csv", dtype={"Data source" : "string", "NCBI tax id":"Int64", "BioSample":"string",
                                          "BioProject":"string","Data Source":"string","NCBI scientific name":"string",
                                          "RDP taxa":"string","RDP taxonomy":"string","basecount":"Int64","16S gene count":"Int64",
                                          "23S gene count":"Int64","tRNA gene count":"Int64","Evidence":"string","Note":"string",
                                          "References":"string"})
df_master = pd.read_csv(path+"/t1d_db_fixed_discussed/good_files_to_use/FINAL_report_WGS_r220_June.csv")

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
species_list_alt = ["Lacrimospora amygdalina","Clostridium sp. Chh4-2","Collinsella aerofaciens","Enterococcus avium","Enterococcus durans",
                    "Lawsonibacter celer","Clostridium symbiosum","Pseudoflavonifractor gallinarum","Roseburia amylophila",
                    "Ruminococcus gnavus","Clostridium perfringens"]

df_master_filtered = df_master[df_master['Selected for Downstream']=="Yes"]
species_list_master = df_master_filtered['GTDBtk Species Classification'].unique().tolist()
species_overlap = []
species_alt_overlap = []

i=0
#print("There are "+str(len(species_list_master))+" unique species in my original list, but the ones missing are: \n")
for item in species_list_master:
    if item not in df['NCBI scientific name'].tolist():
        i=i+1
    else:
        species_overlap.append(item)

i=0
for item in species_list_alt:
    if item not in df['NCBI scientific name'].tolist():
        i=i+1
        #print(item)
    else:
        species_alt_overlap.append(item)
#print("There are "+str(i)+"/"+str(len(species_list_alt))+" species with alternative names missing from rrnDB in the NCBI list:")

list_of_dicts=[]
for list in [species_overlap,species_alt_overlap]:
    overlap_dict={}
    for item in list:
        df_filtered = df[df['NCBI scientific name']==item]
        data_16S = df_filtered['16S gene count']
        median_16S = data_16S.median()
        min_16S = data_16S.min()
        max_16S = data_16S.max()
        count_16S = data_16S.count()
        overlap_dict[item] = {'Count':count_16S,'Median 16S':median_16S,'Minimum 16S':min_16S,'Maximum 16S':max_16S}
    df_16S = pd.DataFrame.from_dict(overlap_dict, orient='index').reset_index()
    list_of_dicts.append(df_16S)

df_rrndb = pd.concat(list_of_dicts).reset_index(drop=True)



#df1_high_q = df1[(df1['CheckM2 Contamination']<5.0) & (df1['CheckM2 Completeness']>95.0) & (df1['High Quality']=="Yes")]
#df1_subset = df1_high_q[['Genome Name / Sample Name', 'IMG Genome ID', 'NCBI Species', 'GTDB Species', 'Assembly Method',
#                  'High Quality', '16S rRNA Count  * assembled', 'CheckM2 Completeness', 'CheckM2 Contamination']]
##df1_subset = df1.groupby('NCBI Species')['16S rRNA Count  * assembled'].agg(['count','mean','median','min','max']).reset_index()
#df1_subset_final = df1.merge(df1_subset, on = 'NCBI Species')
#df1_final = df1_subset_final[['IMG Genome ID','Genome Name / Sample Name','NCBI Species','GTDB Species','High Quality','CheckM2 Completeness','CheckM2 Contamination','count','mean','median','min','max']]
