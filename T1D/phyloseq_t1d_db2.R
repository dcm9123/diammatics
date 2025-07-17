#LOADING PACKAGES
loading_packages = function(){
  library("metagMisc")
  library("decontam")
  library("phyloseq")
  library("ggplot2")
  library("normalize")
  library("dplyr")
  library("ALDEx2")
  library("ANCOMBC")
  library("vegan")
  library("microbiomeMarker")
  library("tidyr")
  library("Bios2cor")
  library("zoo")
  library("patchwork")
  library("Maaslin2")
  library("devtools")
  library("pairwiseAdonis")
  library("DESeq2")
  library("gplots")
  library("reshape2")
  library("Biostrings")
}


#Files needed:
# - seqtab_nochimeras_p1.rds (and for the rest of the 4 plates)
# - plate1/p1_vsearch_dada2_merged.tsv
# - plate1/Nreads_plate1.tsv

###INITIALIZING PATH AND FILE LOCATION
path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/t1d_db_fixed_discussed/FemMicro_Daniel/"
setwd(path)                                                                     #Setting the path of where Im working


###MODIFYING RDS OBJECTS BY ADDING WORD 'PLATE' TO EACH SAMPLE\
#Purpose: Writes a modified RDS file by adding the word 'Plate' to each sample, no other modifications are made
#Input: None
#Returns: Nothing
reshaping_rds_objects = function(){
  obj1 = readRDS("plate1.1/seqtab_nochimeras_p1.rds")
  rownames(obj1) = paste0("Plate1_",rownames(obj1))
  saveRDS(obj1,"plate1/seqtab_nochimeras_m_p1.rds")
  
  obj2 = readRDS("plate2.1/seqtab_nochimeras_p2.rds")
  rownames(obj2) = paste0("Plate2_",rownames(obj2))
  saveRDS(obj2,"plate2/seqtab_nochimeras_m_p2.rds")
  
  obj3 = readRDS("plate3.1/seqtab_nochimeras_p3.rds")
  rownames(obj3) = paste0("Plate3_",rownames(obj3))
  saveRDS(obj3,"plate3/seqtab_nochimeras_m_p3.rds")
  
  obj4 = readRDS("plate4.1/seqtab_nochimeras_p4.rds")
  rownames(obj4) = paste0("Plate4_",rownames(obj4))
  saveRDS(obj4,"plate4/seqtab_nochimeras_m_p4.rds")
  
  obj5 = readRDS("plate5.1/seqtab_nochimeras_p5.rds")
  rownames(obj5) = paste0("Plate5_",rownames(obj5))
  saveRDS(obj5,"plate5/seqtab_nochimeras_m_p5.rds")
} #Adds the word plate followed by its number, that way we create unique IDs for each sample

###RETRIEVING NOCHIMERIC OBJECTS AS RDS FILES
#Purpose: Saves the modified RDS object as seqtab_nochim by plate number
#Input: None
#Returns: Nothing
retrieving_nonchimera_obj = function(){
  seqtab_nochim_p1 <<- readRDS(file = "plate1.1/seqtab_nochimeras_m_p1.rds")  #Reading the non-chimeric sequences from FemMicro using the updated db.
  seqtab_nochim_p2 <<- readRDS(file = "plate2.1/seqtab_nochimeras_m_p2.rds")
  seqtab_nochim_p3 <<- readRDS(file = "plate3.1/seqtab_nochimeras_m_p3.rds")
  seqtab_nochim_p4 <<- readRDS(file = "plate4.1/seqtab_nochimeras_m_p4.rds")
  seqtab_nochim_p5 <<- readRDS(file = "plate5.1/seqtab_nochimeras_m_p5.rds")
} #Gets the original RDS files from FemMicro16S

###INITIALIZING METADATA FILE STRUCTURING
#Purpose: Creates a metadata data frame for the samples in each RDS file
#Input: It takes the seqtab_nochim from the previous function
#Returns: The sample data frame structure (samdf) to be added to each ps object (yet to generate)
initializing = function(variable){
  samples.out<-rownames(variable)                                            #Getting the name of the samples
  names_of_files <- sapply(strsplit(samples.out,"_L001"),`[`,1)                   #This divides the sample names by week, sex, and id, doing it so prevents assigning an S2 or S5 to an illumina sample ID 'S2_L001'
  mice_sex<-logical(length(names_of_files))                                       #Creating a variable to store the sex of the mice
  community_type<-logical(length(names_of_files))                                 #Creating a variable to store the community type (IAB+/-)
  community_subtype<-logical(length(names_of_files))                              #Creating a variable to store the community subtype (NS1,NS6,S2,S5)
  timepoint<-logical(length(names_of_files))
  for (i in seq_along(names_of_files)){                                           #Screening each one of them
    if(grepl("_F_",names_of_files[i]) | grepl("_Female_",names_of_files[i])){
      mice_sex[i]<-"Female"
    }
    else if(grepl("_M_", names_of_files[i]) | grepl("_Male_", names_of_files[i])){
      mice_sex[i]<-"Male"
    }
    else{
      mice_sex[i]<-"Other (check!)"
    }
  }
  
  for (i in seq_along(names_of_files)){                                           #Creating a variable also for mice community and subcommunity
    if(grepl("_S2_",names_of_files[i])){
      community_type[i]<-"IAB+"
      community_subtype[i]<-"S2"
    }
    else if (grepl("_S5_",names_of_files[i])){
      community_type[i]<-"IAB+"
      community_subtype[i]<-"S5"
    }
    else if (grepl("_NS1_",names_of_files[i])){
      community_type[i]<-"IAB-"
      community_subtype[i]<-"NS1"
    }
    else if (grepl("_NS6_",names_of_files[i])){
      community_type[i]<-"IAB-"
      community_subtype[i]<-"NS6"
    }
    else{
      community_type[i]<-"other (check!)"
      community_subtype[i]<-"other (check!)"
    }
  }
  for (i in seq_along(names_of_files)){
    if(grepl("week5", names_of_files[i])){
      timepoint[i]="week_5"
    }
    else if(grepl("week6", names_of_files[i])){
      timepoint[i]="week_6"
    }
    else if(grepl("week7", names_of_files[i])){
      timepoint[i]="week_7"
    }
    else if(grepl("week8", names_of_files[i])){
      timepoint[i]="week_8"
    }
    else if(grepl("week9", names_of_files[i])){
      timepoint[i]="week_9"
    }
    else if(grepl("036R_", names_of_files[i]) & !grepl("week18", names_of_files[i])){
      timepoint[i]="week_9"
    }
    else if(grepl("week10", names_of_files[i])){
      timepoint[i]="week_10"
    }
    else if(grepl("week", names_of_files[i])){
      timepoint[i]=">week_10"
    }
    else if(grepl("diabetes", names_of_files[i])){
      timepoint[i]="week_endpoint"
    }
    else{
      timepoint[i]="NA"
    }
  }
  samdf<-data.frame(ID=samples.out, Sex=mice_sex, Community=community_type, 
                    Subcommunity=community_subtype, Timepoint=timepoint)                               #This generates a data frame by combining all the meatadata I ran the for loop on.
  print(samdf)
  return(samdf)
} #Creates the metadata file by consortia, week, sex, plate number and ID

###INITIALIZING TAXONOMY AND PHYLOSEQ OBJECT GENERATION
#Purpose: Creates a ps object by taking the seqtab_nochim RDS object, the tax table from the FemMicro16S output and the samdf data frame object
#Input: It requires the vsearch_dada2_merged object from FemMicro16S, a seqtab_nochim object previously generated, 
#the samdf file created by the previous function, a prefix for naming it (my case I decided to do p1_, p2_ ...p5_), and
#ASV_sequences = TRUE, so that the sequences become the ID and not the taxa names.
#Returns: The new ps object generated
taxonomy_and_ps = function(file_taxa, seqtab_nochim, samdf,plate_string,asv_sequences){
  taxa_table<-read.table(file_taxa, header=TRUE)                             #This reads the taxonomy table from GTDB by the femmicro pipeline
  taxa_table = taxa_table[,c(2,4,5,6,7,8,9,10)]
  plate_number = strsplit(file_taxa,"/")[[1]][1]
  if(asv_sequences==TRUE){
    asv_ids = colnames(seqtab_nochim)
  }
  else{
    asv_ids = paste0("ASV",seq_len(nrow(taxa_table)))
    print("The ASV numbers will be used as IDs")
    
  }
  taxa_table_with_ids = cbind(ASV_ID = asv_ids, taxa_table)
  write.csv(taxa_table_with_ids,file = paste0(plate_number,"/",plate_number,"_taxa_table.csv"))        #This ensures that the ASV_ID are the same between seqtab and taxa_table after taking the 2nd column                                     #Reading it as a matrix instead of a df
  rownames(taxa_table) = taxa_table[,1]
  taxa_table = taxa_table[,-1]
  rownames(samdf)<-samdf$ID                                                       #Phyloseq doesn't like if the sample names do not match the seqtab file, so renaming the rownames to the ID fixes the problem of "sample names do not match"
  taxa_table<-as.matrix(taxa_table)    
  #print(setdiff(colnames(seqtab_nochim), rownames(taxa_table)))
  ps<-phyloseq(otu_table(seqtab_nochim, taxa_are_rows = FALSE),
               sample_data(samdf),tax_table(taxa_table))
  dna<-Biostrings::DNAStringSet(taxa_names(ps))
  names(dna)<-taxa_names(ps)
  ps<-merge_phyloseq(ps,dna)
  #taxa_names(ps)<-paste0("ASV",plate_string,seq(ntaxa(ps)))
  #print(names(dna))
  #print(ps)
  return(ps)
} #Creates a phyloseq object and ensures the IDs are the same between files when generating one 

###COUNTING ASVs PER RUN
#Purpose: It prints the number of ASVs assigned by GTDB and T1D to the genus and species level and its proportion after cleaning to the class-species level 
#Input: It takes the path and file name of the FemMicro output called vsearch_dada2_merged.tsv object
#Returns: Nothing.
counting_asvs = function(file){
  df = read.table(file, header = TRUE)
  total_asvs = nrow(df)
  kept_asvs = sum(!is.na(df$class_final))
  t1d_asvs = sum(!is.na(df$class_final) & df$database_final=="T1D")
  gtdb_asvs = sum(!is.na(df$class_final) & df$database_final=="GTDB")
  species_asvs = sum(!is.na(df$class_final) & !is.na(df$species_final))
  species_t1d = sum(!is.na(df$class_final) & !is.na(df$species_final) & df$database_final=="T1D")
  species_gtdb = sum(!is.na(df$class_final) & !is.na(df$species_final) & df$database_final=="GTDB")
  genus_asvs = sum(!is.na(df$class_final) & !is.na(df$genus_final))
  genus_t1d = sum(!is.na(df$class_final) & !is.na(df$genus_final) & df$database_final=="T1D")
  genus_gtdb = sum(!is.na(df$class_final) & !is.na(df$genus_final) & df$database_final=="GTDB")
  
  print(paste0("Total ASVs = ",total_asvs))
  print(paste0("Kept ASVs after filtering = ",kept_asvs," (",(round(100*(kept_asvs/total_asvs),2))," %)"))
  print(paste0("ASVs classified with T1D = ",t1d_asvs," (",(round(100*(t1d_asvs/kept_asvs),2))," %)"))
  print(paste0("ASVs classified with GTDB = ",gtdb_asvs," (",(round(100*(gtdb_asvs/kept_asvs),2))," %)"))
  print(paste0("ASVs classified at species-level = ",species_asvs," (",(round(100*(species_asvs/kept_asvs),2))," %)"))
  print(paste0("ASVs speciated with T1D = ",species_t1d," (",(round(100*(species_t1d/species_asvs),2))," %)"))
  print(paste0("ASVs speciated with GTDB = ",species_gtdb," (",(round(100*(species_gtdb/species_asvs),2))," %)"))
  print(paste0("ASVs classified at the genus-level = ",genus_asvs," (",(round(100*(genus_asvs/kept_asvs),2))," %)"))
  print(paste0("ASVs with genus from T1D = ",genus_t1d," (",(round(100*(genus_t1d/genus_asvs),2))," %)"))
  print(paste0("ASVs with genus from GTDB = ",genus_gtdb," (",(round(100*(genus_gtdb/genus_asvs),2))," %)"))
}

###MERGING PS OBJECTS
#Purpose: It merges multiple ps objects
#Input: Five ps objects
#Returns: The merged ps object
merging_runs = function(ps_list){
  ps_merged = do.call(merge_phyloseq,ps_list)
  return(ps_merged)
}

###WRITING ASV AND TAXA TABLE FROM MERGED OBJECT
#Purpose: Writes an ASV count table and a Taxa table to be used for more downstream analyses
#Input: It requires a ps object and a prefix for the file to be written
#Returns: Nothing
writing_tables_ps = function(ps_merged,string){
  merged_asv_table = as.data.frame(t(otu_table(ps_merged))) #transposes the otu table, and transforms it into a data frame
  merged_asv_table$ASV_ID = rownames(merged_asv_table) #Makes a column titled 'ASV_ID'
  merged_asv_table = merged_asv_table[,c(ncol(merged_asv_table), 1:(ncol(merged_asv_table)-1))] #Moves the new ASVID at the beginning
  merged_taxa_table = as.data.frame(tax_table(ps_merged))
  merged_taxa_table$joined = apply(merged_taxa_table, 1, function(x) paste(x, collapse=";")) #Merges all taxonomic ranks and separates them into ';' into one new column 
  sequences = as.character(refseq(ps_merged))
  asv_ids = taxa_names(ps_merged)
  df_sequences = data.frame(ASV_ID = asv_ids,Sequence = sequences, stringsAsFactors = FALSE)
  
  merged_asv_table = merge(merged_asv_table, merged_taxa_table["joined"], by.x="ASV_ID", by.y="row.names")
  merged_asv_table = merge(merged_asv_table, df_sequences, by = "ASV_ID")
  
  merged_taxa_table$Sequence = merged_asv_table$Sequence
  merged_taxa_table$ASV_ID = rownames(merged_taxa_table)
  
  merged_taxa_table <- merged_taxa_table[, c(ncol(merged_taxa_table), 1:(ncol(merged_taxa_table)-1))]  # Move ASV_ID to the front
  write.table(merged_asv_table,paste0(string,"_asv_table.tsv"),sep="\t",row.names = FALSE)
  write.table(merged_taxa_table,paste0(string,"_taxa_table.tsv"), sep="\t",row.names = FALSE)
}

###SUBSETTING MERGED PS BY VARIABLE AND SAMPLE/CTRL TYPE
#Purpose: Creates 4 ps objects. One for the mice samples, positive ctrl, inocula, and negative ctrl. It also counts the total ASVs per sample per subset.
#Input: Takes a merged ps object. In this case, I chose ps_merged_glom (at the taxa level). It will not work unless the actual variable name is passed to the function
#Returns: A list containing the 4 ps objects without any ASVs with a count of zero
subsetting_merged_plate = function(ps_merged_glom){ #I am not sure why the 'ps_object' is not read properly, and I need to pass the 'ps_merged_glom' variable
  ps_merged_weeks = subset_samples(ps_merged_glom, grepl("week",sample_data(ps_merged_glom)$Timepoint, ignore.case=TRUE))
  ps_merged_weeks = prune_taxa(taxa_sums(ps_merged_weeks)>0,ps_merged_weeks)
  weeks_sample_total = sample_sums(ps_merged_weeks)
  write.csv(weeks_sample_total, "sample_total_count_weeks.csv")
  
  ps_merged_positive = subset_samples(ps_merged_glom, grepl("positive|postive", sample_data(ps_merged_glom)$ID, ignore.case=TRUE))
  ps_merged_positive = prune_taxa(taxa_sums(ps_merged_positive)>0,ps_merged_positive)
  positive_sample_total = sample_sums(ps_merged_positive)
  write.csv(positive_sample_total,"sample_total_count_positives.csv")
  
  ps_merged_inocula = subset_samples(ps_merged_glom, grepl("inocul", sample_data(ps_merged_glom)$ID, ignore.case=TRUE))
  ps_merged_inocula = prune_taxa(taxa_sums(ps_merged_inocula)>0, ps_merged_inocula)
  inocula_sample_total = sample_sums(ps_merged_inocula)
  write.csv(inocula_sample_total,"inocula_total_count.csv")
  
  ps_merged_negative = subset_samples(ps_merged_glom, grepl("negative|neg", sample_data(ps_merged_glom)$ID, ignore.case=TRUE))
  ps_merged_negative = prune_taxa(taxa_sums(ps_merged_negative)>0,ps_merged_negative)
  negative_sample_total = sample_sums(ps_merged_negative)
  write.csv(negative_sample_total,"sample_total_count_negatives.csv")
  
  ps_list = list(ps_merged_weeks,ps_merged_inocula,ps_merged_positive,ps_merged_negative)
  return(ps_list)
}

###SUBSETTING WEEKS AND CONSORTIA BY TYPE
#Purpose: To create individual ps objects by week and consortia
#Input: It takes two lists of ps objects. The first one being the ps_merged_weeks, and another one merged by inocula
#Returns: Returns two lists of ps objects divided by week and consortia type
subsetting_weeks_and_consortia = function(ps_weeks){
  ps_by_weeks = list()
  ps_by_consortia = list()
  consortia_names = c("NS1", "NS6", "S2", "S5")
  week_names = c("week_5", "week_6", "week_7", "week_9", "week_10")
  for (i in seq_along(week_names)) {
    week_label = week_names[i]
    idx = sample_data(ps_weeks)$Timepoint == week_label
    ps_by_weeks[[week_label]] = prune_samples(idx, ps_weeks)
    ps_by_weeks[[week_label]] = removing_empty_asvs(ps_by_weeks[[week_label]])
  }
  for (j in seq_along(consortia_names)) {
    cons_label = consortia_names[j]
    idx = sample_data(ps_weeks)$Subcommunity == cons_label
    ps_by_consortia[[cons_label]] = prune_samples(idx, ps_weeks)
    ps_by_consortia[[cons_label]] = removing_empty_asvs(ps_by_consortia[[cons_label]])
  }
  return(list(weeks = ps_by_weeks, consortia = ps_by_consortia))
}

###REMOVING EMPTY ASVs
#Purpose: Removes taxa with no ASV counts
#Input: A ps object
#Returns: A filtered ps object
removing_empty_asvs = function(ps){
  ps_filtered = prune_taxa(taxa_sums(ps)>0, ps)
  return(ps_filtered)
}

###IDENTIFYING OVERLAPPING SAMPLES BY WEEK
#Purpose: Prints the overlapping sample IDs between weeks
#Input: It takes a list of ps objects separated by weeks
#Returns: Nothing
finding_overlap = function(ps_weeks,sample_type){
  duplicates = list()
  combined=list()
  if(sample_type=="weeks"){
    for(i in seq_along(ps_weeks)){
      item = ps_weeks[[i]]
      x = subset_samples(item,Timepoint=="week_5"|Timepoint=="week_6"|Timepoint=="week_7"|Timepoint=="week_9"|Timepoint=="week_10") #Separating ps objects by week
      x = sample_data(x)$ID
      x = sapply(strsplit(x,"_"), function(y) paste(y[2:4], collapse="_"))
      combined[[i]] = x
    }
    matches = list()
    pair_labels = c()
    k = 1  # index for matches list
    
    for (i in 1:(length(combined) - 1)) {
      for (j in (i + 1):length(combined)) {
        intersection = combined[[i]][ combined[[i]] %in% combined[[j]] ]
        matches[[k]] = intersection
        pair_labels[k] = paste0("List", i, "_vs_List", j)
        k = k + 1
      }
    }
  }
  else{
    for(i in seq_along(ps_weeks)){
      item = ps_weeks[[i]]
      x = subset_samples(item,Subcommunity=="NS1"|Subcommunity=="NS6"|Subcommunity=="S2"|Subcommunity=="S5") #Separating ps objects by consortia
      x = sample_data(x)$ID
      x = sapply(strsplit(x,"_"), function(y) paste(y[2:4], collapse="_"))
      combined[[i]] = x
    }
    matches = list()
    pair_labels = c()
    k = 1  # index for matches list
    
    for (i in 1:(length(combined) - 1)) {
      for (j in (i + 1):length(combined)) {
        intersection = combined[[i]][ combined[[i]] %in% combined[[j]] ]
        matches[[k]] = intersection
        pair_labels[k]= paste0("List", i, "_vs_List", j)
        k = k + 1
      }
    }
    for (i in 1:length(combined)){
      duplicates[[i]] = unique(combined[[i]][duplicated(combined[[i]])])
    }
  }
  names(matches) = pair_labels
  print(duplicates)
  return(matches)
}

###PRINT PS INFORMATION
#Purpose: It prints the information pertinent to a ps object
#Input: A ps object
#Returns: Nothing
printing_ps = function(ps){
  print(ps)
}

###CLUSTERING ASVs BY TAXONOMY
#Purpose: To cluster the ASVs to the same taxa if they are pointing like that without removing the NAs
#Input: A ps object and the taxonomic rank to use
#Returns: The clustered ps object
clustering_asvs = function(ps_object,rank){
  ps_clustered = tax_glom(ps_object, taxrank = paste0(rank,"_final"), NArm = FALSE)
  return(ps_clustered)
}

###RENAMING THE ASV IDs (in this case sequences) INTO ITS TAXONOMY
#Purpose: Renames the ASV IDs into the actual taxa that was previously clustered
#Input: A tax_glom ps object
#Return: The modified ps object
otu_ID_renaming = function(ps_object){
  new_name = paste0(tax_table(ps_object)[,"genus_final"],"_",tax_table(ps_object)[,"species_final"])
  new_name = as.character(new_name)
  new_name = make.unique(new_name)
  taxa_names(ps_object) = new_name
  length
  return(ps_object)
}

###FILTERING PS OBJECT BY SPECIFIC SAMPLES
#Purpose: To remove specific samples passed by a vector
#Input: A vector of sample IDs and the list of ps objects they will be removed from
#Output: A ps object without those samples
filtering_samples = function(id_vector, ps_obj2){
  chopped_IDs = sapply(strsplit(as.character(sample_data(ps_obj2)$ID),"_"),function(x) paste(x[2:4], collapse="_"))
  cage_id = logical(length(chopped_IDs))
  discarding = chopped_IDs %in% id_vector
  keeping = sample_names(ps_obj2)[!(discarding)]
  print(paste0("The original number of samples in the ps object: ", nsamples(ps_obj2)))
  ps_obj2 = prune_samples(keeping, ps_obj2)
  print(paste0("The kept number of samples in the ps object: ", nsamples(ps_obj2)))
  return(ps_obj2)
}

#WRITING AN OTU/COUNT TABLE
#Purpose: Writing a otu table in .csv format from the ps object that is passed
#Input: A ps object
#Return: none
otu_writing = function(ps_obj,name_out){
  otu_t = otu_table(ps_obj)
  write.csv(otu_t,file = paste0(name_out))
}

#################################################################################
#                                                                               #
#                               MOTHER BOARD                                    #
#                                                                               #
#                                                                               #
#################################################################################
path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/t1d_db_fixed_discussed/FemMicro_Daniel/"
setwd(path)  
loading_packages()
reshaping_rds_objects()
retrieving_nonchimera_obj() #This function defines the seqtab_nochim variables inside of it and sends it to the global env.

#Getting samdfs for each ps object (yet to be built)
samdf1 = initializing(seqtab_nochim_p1)
samdf2 = initializing(seqtab_nochim_p2)
samdf3 = initializing(seqtab_nochim_p3)
samdf4 = initializing(seqtab_nochim_p4)
samdf5 = initializing(seqtab_nochim_p5)

#Creating ps objects, leaving asv_sequences = TRUE to get ASV sequences as IDs and not taxa names
ps1 = taxonomy_and_ps("plate1.1/final_merged_tables/p1_vsearch_dada2_merged.tsv",seqtab_nochim_p1,samdf1,"p1_",asv_sequences = TRUE)
ps2 = taxonomy_and_ps("plate2.1/final_merged_tables/p2_vsearch_dada2_merged.tsv",seqtab_nochim_p2,samdf2,"p2_",asv_sequences = TRUE)
ps3 = taxonomy_and_ps("plate3.1/final_merged_tables/p3_vsearch_dada2_merged.tsv",seqtab_nochim_p3,samdf3,"p3_",asv_sequences = TRUE)
ps4 = taxonomy_and_ps("plate4.1/final_merged_tables/p4_vsearch_dada2_merged.tsv",seqtab_nochim_p4,samdf4,"p4_",asv_sequences = TRUE)
ps5 = taxonomy_and_ps("plate5.1/final_merged_tables/p5_vsearch_dada2_merged.tsv",seqtab_nochim_p5,samdf5,"p5_",asv_sequences = TRUE)

#Counting reads by plate
#counting_reads("plate1.1/Nreads_plate1.tsv")
#counting_reads("plate2.1/Nreads_plate2.tsv")
#counting_reads("plate3.1/Nreads_plate3.tsv")
#counting_reads("plate4.1/Nreads_plate4.tsv")
#counting_reads("plate5.1/Nreads_plate5.tsv")

#Counting ASVs by plate
counting_asvs("plate1.1/final_merged_tables/p1_vsearch_dada2_merged.tsv")
counting_asvs("plate2.1/final_merged_tables/p2_vsearch_dada2_merged.tsv")
counting_asvs("plate3.1/final_merged_tables/p3_vsearch_dada2_merged.tsv")
counting_asvs("plate4.1/final_merged_tables/p4_vsearch_dada2_merged.tsv")
counting_asvs("plate5.1/final_merged_tables/p5_vsearch_dada2_merged.tsv")

#Printing raw ps object generated
printing_ps(ps1) #172 ASVs and 38 samples
printing_ps(ps2) #161 ASVs and 90 samples
printing_ps(ps3) #187 ASVs and 90 samples
printing_ps(ps4) #279 ASVs and 90 samples
printing_ps(ps5) #102 ASVs and 56 samples

#Merging raw ps objects into one
ps_list1 = list(ps1,ps2,ps3,ps4,ps5)
ps_merged = merging_runs(ps_list1)
printing_ps(ps_merged) #535 ASVs and 364 samples (535 after removing identical sequences, 901 if counted separately by ps object)

#Clustering ASVs to species level and not removing the NAs
ps_clustered = clustering_asvs(ps_merged,"species")
printing_ps(ps_clustered) #99 species (or taxa) and 364 samples

#Subsetting the merged plate into sample weeks, inocula, positives, and negatives and renaming the OTU IDs by taxonomy, not by sequence
ps_clustered_list = subsetting_merged_plate(ps_clustered) #1: weeks, #2: inocula, #3: positive, #4: negative
ps_clustered_renamed = list()
for(i in 1:4){
  ps_clustered_renamed[[i]] = otu_ID_renaming(ps_clustered_list[[i]]) #Calling the renaming function
  printing_ps(ps_clustered_renamed[[i]])
}

#Dividing ps object by week and consortia
ps_weeks = ps_clustered_renamed[[1]] #ps dividing by mice samples only
ps_weeks_and_consortia_double_list = subsetting_weeks_and_consortia(ps_weeks)

#Assigning cleaned ps corresponsing variable
ps_w5 = ps_weeks_and_consortia_double_list$weeks$week_5 #week 5 55 species and 115 samples
ps_w6 = ps_weeks_and_consortia_double_list$weeks$week_6 #week 6 32 species and 16 samples
ps_w7 = ps_weeks_and_consortia_double_list$weeks$week_7 #week 7 31 species and 16 samples
ps_w9 = ps_weeks_and_consortia_double_list$weeks$week_9 #week 9 54 species and 119 samples
ps_w10 = ps_weeks_and_consortia_double_list$weeks$week_10 #week 10 39 species and 29 samples
ps_ns1 = ps_weeks_and_consortia_double_list$consortia$NS1 #ns1 39 species and 127 samples
ps_ns6 = ps_weeks_and_consortia_double_list$consortia$NS6 #ns6 31 species and 32 samples
ps_s2 = ps_weeks_and_consortia_double_list$consortia$S2 #s2 46 species and 148 samples
ps_s5 = ps_weeks_and_consortia_double_list$consortia$S5 #s5 38 species and 32 samples

#There are a total of 295 mice with week identifiers in them, the rest belong into the cage and diabetes endpoint analysis
for(item in ps_weeks_and_consortia_double_list){
  printing_ps(item)
}

#This section here will identify the repeated samples:
overlapping_samples_by_weeks = finding_overlap(ps_weeks_and_consortia_double_list[[1]], "weeks") #[[1]] is the mice samples divided by week
overlapping_samples_by_consortia = finding_overlap(ps_weeks_and_consortia_double_list[[2]], "consortia") #[[2]] is the mice samples divided by consortia

overlapping_samples_by_weeks

for(item in sample_data(ps_s2)$ID){ #16 Samples, all week 6 are present in week 5, but no repetitions are found between week 9 and week 10
  cat(item)
  cat("\n")
}

#Filtering samples by repetitions in consortia
ps_w5_f = filtering_samples(overlapping_samples_by_weeks[[1]],ps_w5) #This command removes the duplicates of the same mice across w5 and w6
                                                                     #There was no point in running the rest, as there were no other repetitions found
ps_w5_f
ps_list2 = list(ps_w5_f,ps_w6,ps_w7,ps_w9,ps_w10)
ps_merged2 = merging_runs(ps_list2) #This has a total of 279 samples, which is the original 295 - 16 repetitions from w6 in w5
ps_merged2_divided = subsetting_weeks_and_consortia(ps_merged2) #This divides the ps into weeks [[1]] and consortia [[2]]
ps_consortia_filtered = ps_merged2_divided[[2]] #109 to 101 from NS1, 122 to 114 in S2 reduction, the 16 repeated samples are here
ps_ns1_final = ps_consortia_filtered[[1]]
ps_ns6_final = ps_consortia_filtered[[2]]
ps_s2_final = ps_consortia_filtered[[3]]
ps_s5_final = ps_consortia_filtered[[4]]

#Writing out the otu table for week 5, week 9, and week 10
write.csv(t(otu_table(ps_ns1_final)),"ps_ns1_final.csv")
write.csv(t(otu_table(ps_ns6_final)),"ps_ns6_final.csv")
write.csv(t(otu_table(ps_s2_final)),"ps_s2_final.csv")
write.csv(t(otu_table(ps_s5_final)),"ps_s5_final.csv")

