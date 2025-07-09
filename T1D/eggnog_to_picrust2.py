# Daniel Castaneda Mogollon, PhD
# 14:49
# eggnog_to_picrust.py.py
# Purpose: This script was made to format the output from emapper to something picrust2 can understand.
#It is inteded to format EC and KOs.

import pandas as pd
import os
import sys
from tqdm import tqdm


def setting_up():
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 0)  # Let pandas auto-detect width

def welcome():
    #This function welcomes the user and asks whether they should run the code for an individual folder (with no merging) or multiple folders and a final merged file.
    print("This code has been created to transform the format from an eggnog_mapper.annotation file to be readable by"
          "PICRUSt2. It is IMPERATIVE that the name of the samples are in this format: S_NS1_Blss_094.emapper.annotations."
          "First instance is S identifier, second is the consortia it belongs to, third it's the original species abbreviation,"
          "and fourth, is the sample number. You only need to change the directory where these files are located in the 'path' "
          "variable at the end.")
    list_paths=[]
    response = input("Would you like to run the entire code across MULTIPLE folders to merge into a final file? Or just an individual folder to format the eggnog file to PICRUSt2? (yes/no)")
    if response.lower() == "no":
        answer = input("Please provide the path where your eggnog files are located: ")
        list_paths.append(answer)
    elif response.lower() == "yes":
        number_of_paths = input("Please provide the number of paths you want to work with: ")
        number_of_paths = int(number_of_paths)
        list_paths=[]
        for number in range(0,number_of_paths,1): #This iterates over the paths the user will give
            path = input("Please provide your #" +str(number+1)+" of where your eggnog files are located: ")
            if not os.path.exists(path):
                print("Path doesn't exist. Terminating program now")
                sys.exit()
            list_paths.append(path)
    else:
        print("Not a valid answer. Terminating program now.")
        sys.exit()
    answer_output = input("Please provide the OUTPUT path where you want your eggnog formatted files to be written to "
                          "(if it doesn't exist, we'll create one for you):")
    os.makedirs(answer_output, exist_ok=True) #Creating a dir in case it isn't there
    return(list_paths,answer_output) #This returns the input path(s) and the output path

#Calling the first function to format the original eggnog file to something pandas-readable
def formatting_original_file(path,answer_output): #This function formats the eggnog file by removing unnecessary lines and adding columns
    os.chdir(path)
    consortia_names = ["NS1","NS6","S2","S5"]

    total_files = os.listdir()
    pbar = tqdm(total=len(total_files),desc="Reading files from "+path)

    for folder in consortia_names:
        os.makedirs(answer_output+"/"+folder, exist_ok=True)
    for file in os.listdir(path):
        if file.endswith(".annotations"):
            prefix = file.split("_")[1]  # Removes the _eggnog.tsv as part of the variable. This will be used to name the genomes
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
                df.to_csv(answer_output+'/'+prefix+'/'+file_name + "_eggnog.tsv", sep='\t', index=False)  # Output as tsv file
                #print(f"{file_name} has been written to {answer_output}/{prefix}_{file_name}_eggnog.tsv")
            os.remove(f_out)
        else:
            print(file+" is not an emapper.annotations file, skipping...")
        pbar.update(1)
    pbar.close()
    return(answer_output)

def formatting_picrust2_annotations(path):
    #This function is the meat of the code. It takes the ECs and KOs out of the previously generated file
    #to extract the ECs and KOs, count them (also removes duplicates), and write them in a new .tsv file for PICRUSt2
    #it also reformats the way the KOs in KEGG are writen, and the ECs too.
    os.chdir(path)
    consortia_names = ["NS1","NS6","S2","S5"]
    total_files = 0
    for consortium in consortia_names:
        for file in os.listdir(path+"/"+consortium):
            if file.endswith("_eggnog.tsv"):
                total_files = total_files + 1

    pbar = tqdm(total=total_files, desc="Formatting files for PICRUSt2")    #For loading bar

    for consortium in consortia_names:
        os.chdir(path+'/'+consortium)
        for file in os.listdir():
            if file.endswith("_eggnog.tsv"): #Works with the newly generated files from the previous function.
                prefix = file.split("_eggnog")[0] #Removes the _eggnog.tsv as part of the variable. This will be used to name the genomes
                df = pd.read_csv(file,sep='\t') #Reading the file
                df["EC"] = df["EC"].replace("-","NA")
                df["EC"] = df["EC"].astype(str).str.replace(r'(\d+(?:\.\d+)+)', r'EC:\1', regex=True)
                df["KEGG_ko"] = df["KEGG_ko"].replace("-","NA") #Replaces entire cell containing '-'
                df["KEGG_ko"] = df["KEGG_ko"].str.replace("ko:","") #Replaces instances of "ko:" in a cell, that's why we used str.replace
                df_subset = df.iloc[:,[0,11,12]].copy() #This calls only the ID of the sample, the ECs, and the KOs, and makes sure I work on a copy and do not modify the original df
                for column in df_subset.columns[1:3]:
                    os.makedirs(path+"/"+consortium+"/"+column.replace("KEGG_",""), exist_ok=True)
                    df_subset[column] = df_subset[column].astype(str) #Makes sure that the EC numbers are treated as a string
                    df_subset[column] = df_subset[column].str.split(",") #Separates the multiple EC numbers of each cell into different rows
                    df_exploded_ec_ko = df_subset.explode(column) #Separates into new rows the split values
                    ec_ko_counts = df_exploded_ec_ko[column].value_counts() #Counts the instances of each EC value
                    ec_ko_new = pd.DataFrame([ec_ko_counts.to_dict()]) #Turns these into a dictionary, where the key is the EC number, and the value is the count. It is transformed from series to df
                    ec_ko_new.drop(columns="NA", inplace=True) #Deletes the 'NA' value column
                    ec_ko_new.insert(0, "Sample ID",prefix)
                    #print("Writing file. . . "+path+"/"+consortium+"/"+prefix+"_"+column.replace("KEGG_","")+"_eggnog.tsv")
                    ec_ko_new.to_csv(path+"/"+consortium+"/"+str(column.replace("KEGG_",""))+"/"+prefix+"_"+column.replace("KEGG_","")+"_eggnog.tsv", sep='\t', index=False) #Writes the sample in its corresponding folder.
                pbar.update(1)
    pbar.close()

def merging(path):
    df_list_ko=[]
    df_list_ec=[]
    consortia = ["NS1","NS6","S2","S5"]
    annotations = ["ko","EC"]

    total_files = 0 #This section was created to implement a loading bar
    for consortium in consortia:
        for annotation in annotations:
            total_files = total_files+len(os.listdir(path+"/"+consortium+"/"+annotation)) #This gives me the final number of files, so I can use it for the loading bar

    pbar = tqdm(total=total_files, desc="Merging and writing all the files for KO and EC annotations")

    for consortium in consortia:
        for annotation in annotations:
            new_path = path+"/"+consortium+"/"+annotation
            os.chdir(new_path)
            for file in os.listdir():
                if annotation=="ko":
                    df_list_ko.append(pd.read_csv(new_path+'/'+file,sep='\t'))
                else:
                    df_list_ec.append(pd.read_csv(new_path+'/'+file,sep='\t'))
                pbar.update(1)

    df_final_ko = pd.concat(df_list_ko,axis=0,join="outer")
    df_final_ec = pd.concat(df_list_ec,axis=0,join="outer")
    df_final_ko.to_csv(path+"/Final_KO_file.tsv",sep='\t',index=False)
    df_final_ec.to_csv(path+"/Final_EC_file.tsv",sep='\t',index=False)
    pbar.close()


def main():
    list_of_paths,output_path = welcome() #We get a tupple here, so I call these two separately
    for path in list_of_paths:
        formatting_original_file(path,output_path)
    formatting_picrust2_annotations(output_path)
    merging(output_path)
    print("Done!")
if __name__ == "__main__":
    main()
