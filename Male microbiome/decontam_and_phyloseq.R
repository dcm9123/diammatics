#############MICROBIOME ANALYSIS#########################
#                                                       #
#                                                       #
#########################################################

#Author: Daniel Castaneda Mogollon, PhD
#Date: March 22nd, 2024
#Sycuro Lab

###INSTALLING PACKAGES###

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("microbiomeMarker")
BiocManager::install("decontam")
BiocManager::install("normalize")
BiocManager::install("ALDEx2")
BiocManager::install("ANCOMBC")
BiocManager::install("betatest")
remotes::install_github("vmikk/metagMisc")
BiocManager::install("Bios2cor")
BiocManager::install("Maaslin2")
BiocManager::install("picante")
install.packages("RAxML")

#LOADING PACKAGES

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
library("picante")
library("raxml")

#IMPORTING OBJECTS AND SETTING DIRECTORIES

path = "/Users/danielcm/Desktop/Sycuro/Projects/Male_microbiome/male_microbiome/v1v2/output/dada2"
setwd(path)
seqtab.nochim <- readRDS(file = "seqtab_nochimeras.rds")
samples.out <- rownames(seqtab.nochim)
samples.out
names_of_files <- sapply(strsplit(samples.out,"_"),`[`,1) #This divides the sample sites as whole, penile, urine, Undetermined, ATCC
contains_dummy <- logical(length(names_of_files))
for (i in seq_along(names_of_files)){
  if (grepl("Dummy", names_of_files[i])){
    contains_dummy[i]<-"Dummy"
  }
  else if(grepl("ATCC", names_of_files[i]) | grepl("Amplification",names_of_files[i])){
    contains_dummy[i]<-"Control"
  }
  else{
    contains_dummy[i]<-"Sample"
  }
}
print(contains_dummy)
subset_names_of_files<-sapply(strsplit(names_of_files,"-"),`[`,1)
samples_id <-paste(subset_names_of_files,contains_dummy,sep="-")
kit_type<-logical(length(names_of_files))  #This will subset each file by its kit type 
for (i in seq(names_of_files)){
  if (grepl("Norgen", names_of_files[i])){
    kit_type[i]<-"Norgen"
  }
  else if (grepl("Bact", names_of_files[i])){
    kit_type[i]<-"Bact"
  }
  else if (grepl("PSP", names_of_files[i])){
    kit_type[i]<-"PSP"
  }
  else if (grepl("Skin", names_of_files[i])){
    kit_type[i]<-"PSP"
  }
  else{
    kit_type[i]<-"Unknown"
  }
}

sample_type<-logical(length(names_of_files))
names_of_files
for (i in seq(names_of_files)){
  if (grepl("Skin", names_of_files[i])){
    sample_type[i]<-"Penile-Skin"
  }
  else if (grepl("Pellet", names_of_files[i])){
    sample_type[i]<-"Urine-Pellet"
  }
  else if (grepl("Whole-Urine", names_of_files[i])){
    sample_type[i]<-"Whole-Urine"
  }
  else{
    sample_type[i]<-"Other"
  }
}
sample_type
samdf <-data.frame(ID=samples.out, Type=samples_id, Type2=sample_type, Subtype=contains_dummy, Kit=kit_type)
samdf
rownames(samdf)<-samples.out


#COUNTING RAW AND FILTERED READS FROM DADA2

reads_table<-as.data.frame(read.table("Nreads.tsv"))
total_raw_reads<-sum(as.numeric(reads_table$V2[-1]))
total_nochim_reads<-sum(as.numeric(reads_table$V7[-1]))
psp_skin_raw_reads<-sum(as.numeric(reads_table$V2[2:31]))
psp_skin_filtered_reads<-sum(as.numeric(reads_table$V7[2:31]))
bact_raw_reads<-sum(as.numeric(reads_table$V2[33:38]))
bact_filtered_reads<-sum(as.numeric(reads_table$V7[33:38]))
psp_urine_raw_reads<-sum(as.numeric(reads_table$V2[39:65]))
psp_urine_filtered_reads<-sum(as.numeric(reads_table$V7[39:65]))
norgen_raw_reads<-sum(as.numeric(reads_table$V2[66:98]))
norgen_filtered_reads<-sum(as.numeric(reads_table$V7[66:98]))
reads_table
reads_table$V1[2:31]

print(total_raw_reads)
print(total_nochim_reads)
print(norgen_raw_reads)
print(norgen_filtered_reads)
print(psp_skin_raw_reads)
print(psp_skin_filtered_reads)
print(psp_urine_raw_reads)
print(psp_urine_filtered_reads)
print(sum(psp_skin_raw_reads,psp_urine_raw_reads))
print(sum(psp_skin_filtered_reads,psp_urine_filtered_reads))
print(bact_raw_reads)
print(bact_filtered_reads)


###TAXA TABLE
taxa_table<-read.table("../taxonomy/GTDB_RDP.tsv",header=TRUE)
taxa_table<-as.matrix(taxa_table) #Sometimes this works and sometimes it doesnt...
#View(taxa_table)
taxa_table

###PHYLOSEQ OBJECT
ps<-phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
                                  sample_data(samdf),
                                  tax_table(taxa_table))
dna<-Biostrings::DNAStringSet(taxa_names(ps))
names(dna)<-taxa_names(ps)
ps<-merge_phyloseq(ps,dna)
taxa_names(ps)<-paste0("ASV",seq(ntaxa(ps)))

otu_table(ps)
sample_data(ps)
tax_table(ps)

#COUNTING THE NUMBER OF TAXA IN THE PS OBJECT
ps
sum(colSums(ps@otu_table)>1)
length(taxa_names(ps))

###SUBSETTING EACH PS OBJECT BY KIT. 
#Each Kit was a different run, and therefore should be decontaminated separately
ps_norgen<-subset_samples(ps, Kit=="Norgen")
ps_psp<-subset_samples(ps, Kit=="PSP")
ps_bact<-subset_samples(ps, Kit=="Bact")

sample_data(ps_norgen)
sample_data(ps_psp)

#PRUNNING TAXA IF THEY ARE ABSCENT BY KIT
ps_norgen_pruned<-prune_taxa(taxa_sums(ps_norgen)>0, ps_norgen)
ps_psp_pruned<-prune_taxa(taxa_sums(ps_psp)>0, ps_psp)
ps_bact_pruned<-prune_taxa(taxa_sums(ps_bact)>0, ps_bact)
ps_norgen_pruned
ps_psp_pruned
ps_bact_pruned
sample_data(ps_psp_pruned)
#Now dividing by the body site of the same kit, psp and bact (norgen to be ommited... for now)
ps_psp_skin_pruned<-prune_samples(sample_data(ps_psp_pruned)$Type2=="Penile-Skin", ps_psp_pruned)
ps_psp_urine_pruned<-prune_samples(sample_data(ps_psp_pruned)$Type2=="Urine-Pellet", ps_psp_pruned)
ps_psp_skin_pruned<-prune_taxa(taxa_sums(ps_psp_skin_pruned)>0, ps_psp_skin_pruned)
ps_psp_urine_pruned<-prune_taxa(taxa_sums(ps_psp_urine_pruned)>0, ps_psp_urine_pruned)

ps_psp_skin_pruned
ps_psp_urine_pruned

#Leaving only the samples in a ps object (removing controls and undetermined)
ps_psp_skin_samples_only<-prune_samples(sample_data(ps_psp_skin_pruned)$Type=="Penile-Sample", ps_psp_skin_pruned)
ps_psp_skin_samples_only<-prune_taxa(taxa_sums(ps_psp_skin_samples_only)>0, ps_psp_skin_samples_only)
ps_psp_urine_samples_only<-prune_samples(sample_data(ps_psp_urine_pruned)$Type=="Urine-Sample", ps_psp_urine_pruned)
ps_psp_urine_samples_only<-prune_taxa(taxa_sums(ps_psp_urine_samples_only)>0, ps_psp_urine_samples_only)
ps_bact_samples_only<-prune_samples(sample_data(ps_bact_pruned)$Type=="Urine-Sample", ps_bact_pruned)
ps_bact_samples_only<-prune_taxa(taxa_sums(ps_bact_samples_only)>0, ps_bact_samples_only)

ps_bact_samples_only
ps_psp_skin_samples_only
ps_psp_urine_samples_only

###RUNNING POSTIIVE CONTROL ATCC 
sample_data(ps)
ps_atcc<-prune_samples(sample_data(ps)$Subtype=="Control", ps)
ps_atcc<-prune_taxa(taxa_sums(ps_atcc)>0, ps_atcc)
otu_table(ps_atcc)
tax_table(ps_atcc)


####DECONTAM PACKAGE BELOW
df_sampledata<-as.data.frame(sample_data(ps))
df_sampledata$librarysize<-sample_sums(ps)
df_sampledata<-df_sampledata[order(df_sampledata$librarysize)]
df_sampledata$index<-seq(nrow(df_sampledata))
for (value in df_sampledata$librarysize){
  cat(value,"\n")
}
ggplot(data=df_sampledata, aes(x=index, y=librarysize, color=Subtype)) + geom_point()
df_sampledata


#PREVALENCE METHOD FOR DECONTAMINATION (FREQUENCY IS WHEN WE HAVE THE DNA CONCENTRATION OF ALL SAMPLES VIA A FLUORESCENCE ASSAY)
#REPEAT THIS STEP FOR EACH KIT TYPE (i.e. ps_psp, ps_norgen, ps_bact) MAKING NEW PS OBJECTS WITHOUT CONTAMINANTS
thres = 0.5 #User-defined. A 0.5 means that all the ASVs found in the control higher than the samples will be removed from the samples themselves.

#Qiagen bact kit first
sample_data(ps_bact)$is.neg<-sample_data(ps_bact)$Subtype == "Dummy"
contamdf.prev<-isContaminant(ps_bact, method="prevalence", neg="is.neg",threshold=thres) #threshold of 0.5 identifies all contaminants as all sequences that are more prevalent in negative control than in the positive samples
table(contamdf.prev$contaminant)
which(contamdf.prev$contaminant) #This prints the ASV number of my contaminant in this data set.
bact_contaminant<-which(contamdf.prev$contaminant)
bact_contaminant
class(bact_contaminant)
if(identical(bact_contaminant,integer(0))){
  print("No ASV founds in this variable ps")
}else{
  bact_contaminant<-paste0("ASV", bact_contaminant)
  tax_table(ps_bact)[bact_contaminant,]
  }

bact_contaminant
ps_bact_clean<-subset_taxa(ps_bact, !taxa_names(ps_bact) %in% bact_contaminant)
ps.pa<-transform_sample_counts(ps_bact, function(abund) 1*(abund>0))
ps.pa.neg<-prune_samples(sample_data(ps.pa)$Subtype=="Dummy",ps.pa)
ps.pa.pos<-prune_samples(sample_data(ps.pa)$Subtype=="Sample",ps.pa)
df.pa<-data.frame(pa.pos=taxa_sums(ps.pa.pos),pa.neg=taxa_sums(ps.pa.neg), contaminant=contamdf.prev$contaminant)
df.pa
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) +geom_point() + xlab("Prevalence (Dummies)") +ylab("Prevalence (Samples)") + ggtitle("v2v3") + theme(plot.title = element_text(hjust = 0.5))
contamdf.prev
sample_data(ps.pa.pos)

bact_decontaminated<-prune_samples(sample_data(ps_bact_clean)$Type=="Urine-Sample",ps_bact_clean)
bact_decontaminated<-prune_taxa(taxa_sums(bact_decontaminated)>0, bact_decontaminated)
bact_decontaminated<-prune_samples(sample_sums(bact_decontaminated)>0, bact_decontaminated)

#Norgen kit next
sample_data(ps_norgen)$is.neg<-sample_data(ps_norgen)$Subtype == "Dummy"
contamdf.prev<-isContaminant(ps_norgen, method="prevalence", neg="is.neg",threshold=thres) #threshold of 0.5 identifies all contaminants as all sequences that are more prevalent in negative control than in the positive samples
table(contamdf.prev$contaminant)
which(contamdf.prev$contaminant) #This prints the ASV number of my contaminant in this data set.
norgen_contaminant<-which(contamdf.prev$contaminant) 
norgen_contaminant<-paste0("ASV",norgen_contaminant) #Assigning 'ASV' to the contaminant number
tax_table(ps_norgen)[norgen_contaminant]
ps_norgen_clean<-subset_taxa(ps_norgen, !taxa_names(ps_norgen) %in% norgen_contaminant) #Removing contaminants from the data set, and making a clean new one
ps_norgen_clean #Sanity check
ps.pa<-transform_sample_counts(ps_norgen, function(abund) 1*(abund>0))
ps.pa.neg<-prune_samples(sample_data(ps.pa)$Subtype=="Dummy",ps.pa)
ps.pa.pos<-prune_samples(sample_data(ps.pa)$Subtype=="Sample",ps.pa)
df.pa<-data.frame(pa.pos=taxa_sums(ps.pa.pos),pa.neg=taxa_sums(ps.pa.neg), contaminant=contamdf.prev$contaminant)
df.pa
ps.pa.pos
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) +geom_point() + xlab("Prevalence (Dummies)") +ylab("Prevalence (Samples)") + ggtitle("v2v3") + theme(plot.title = element_text(hjust = 0.5))

#PSP kit last
sample_data(ps_psp)$is.neg<-sample_data(ps_psp)$Subtype == "Dummy"
contamdf.prev<-isContaminant(ps_psp, method="prevalence", neg="is.neg",threshold=thres) #threshold of 0.5 identifies all contaminants as all sequences that are more prevalent in negative control than in the positive samples
table(contamdf.prev$contaminant)
which(contamdf.prev$contaminant) #This prints the ASV number of my contaminant in this data set.
psp_contaminant<-which(contamdf.prev$contaminant)
psp_contaminant<-paste0("ASV", psp_contaminant)
tax_table(ps_psp)[psp_contaminant]
ps_psp_clean<-subset_taxa(ps_psp, !taxa_names(ps_psp) %in% psp_contaminant)
ps_psp_clean #sanity check
ps.pa<-transform_sample_counts(ps_psp, function(abund) 1*(abund>0))
ps.pa.neg<-prune_samples(sample_data(ps.pa)$Subtype=="Dummy",ps.pa)
ps.pa.pos<-prune_samples(sample_data(ps.pa)$Subtype=="Sample",ps.pa)
df.pa<-data.frame(pa.pos=taxa_sums(ps.pa.pos),pa.neg=taxa_sums(ps.pa.neg), contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) +geom_point() + xlab("Prevalence (Dummies)") +ylab("Prevalence (Samples)") + ggtitle("v2v3") + theme(plot.title = element_text(hjust = 0.5))
contamdf.prev

psp_skin_decontaminated<-prune_samples(sample_data(ps_psp_clean)$Type=="Penile-Sample", ps_psp_clean) #No controls
psp_skin_decontaminated<-prune_taxa(taxa_sums(psp_skin_decontaminated)>0, psp_skin_decontaminated) #No empty taxa included
psp_skin_decontaminated<-prune_samples(sample_sums(psp_skin_decontaminated)>0, psp_skin_decontaminated)

psp_urine_decontaminated<-prune_samples(sample_data(ps_psp_clean)$Type=="Urine-Sample",ps_psp_clean)
psp_urine_decontaminated<-prune_taxa(taxa_sums(psp_urine_decontaminated)>0, psp_urine_decontaminated)
psp_urine_decontaminated<-prune_samples(sample_sums(psp_urine_decontaminated)>0, psp_urine_decontaminated)

#Estimating the final taxa by ps decontaminated object and only including samples
psp_skin_decontaminated
psp_urine_decontaminated
bact_decontaminated


#Writing the PS otu table objects before decontaminating (aka raw)
write.csv(file = "ps_psp.csv", x = otu_table(ps_psp))
write.csv(file = "ps_psp_skin.csv", x = otu_table(ps_psp_skin_samples_only))
write.csv(file = "ps_psp_urine.csv", x = otu_table(ps_psp_urine_samples_only))
write.csv(file = "ps_norgen.csv", x = otu_table(ps_norgen))
write.csv(file = "ps_bact.csv", x = otu_table(ps_bact))

#Writing the decontaminated ASV table files from each kit and sample type
write.csv(file = "ps_psp_skin_clean.csv", x = otu_table(psp_skin_decontaminated))
write.csv(file = "ps_psp_urine_clean.csv", x = otu_table(psp_urine_decontaminated))
write.csv(file = "ps_norgen_clean.csv", x = otu_table(ps_norgen_clean))
write.csv(file = "ps_bact_clean.csv", x = otu_table(bact_decontaminated))

#GETTING FILES TO PLOT BEFORE AND AFTER CONTAMINATION BY KIT
sample_sums_bact<-sample_sums(ps_bact)
sample_sums_bact_clean<-sample_sums(ps_bact_clean)
sample_sums_psp<-sample_sums(ps_psp)
sample_sums_psp_clean<-sample_sums(ps_psp_clean)
sample_sums_norgen<-sample_sums(ps_norgen)
sample_sums_norgen_clean<-sample_sums(ps_norgen_clean)
cat(paste0(sample_sums_bact,"\n"))

cat(unname(paste0(sample_sums_bact,"\n")))
cat(unname(paste0(sample_sums_bact_clean,"\n")))
cat(unname(paste0(sample_sums_psp,"\n")))
cat(unname(paste0(sample_sums_psp_clean,"\n")))
cat(unname(paste0(sample_sums_norgen,"\n")))
cat(unname(paste0(sample_sums_norgen_clean,"\n")))


#ESTIMATING THE NUMBER OF PHYLA, CLASS, ORDER, FAMILY, GENUS, AND SPECIES ACROSS REGION AND SAMPLES
taxa_psp_skin<-as.data.frame(tax_table(psp_skin_decontaminated))
taxa_psp_urine<-as.data.frame(tax_table(psp_urine_decontaminated))
taxa_bact_urine<-as.data.frame(tax_table(bact_decontaminated))
taxa_ps_objects<-c(taxa_psp_skin, taxa_psp_urine, taxa_bact_urine)
taxonomic_ranks<-c("Phylum","Class","Order","Family","Genus","Species")
taxa_psp_urine$Phylum

print("PSP SKIN DATA:")
for(rank in taxonomic_ranks){
  rank_psp_skin<-taxa_psp_skin[, rank][!is.na(taxa_psp_skin)]
  distinct_rank_psp_skin<-length(unique(rank_psp_skin))
  total_taxa_psp_skin<-nrow(taxa_psp_skin)
  length(taxa_psp_skin[, rank][is.na(taxa_psp_skin[, rank])])
  total_taxa_psp_skin
  na_count_psp_skin<-length(taxa_psp_skin[, rank][is.na(taxa_psp_skin[, rank])])
  proportion_na_psp_skin<-na_count_psp_skin/total_taxa_psp_skin
  print(paste("The number of distinct",rank,"(excluding NAs) in the PSP skin samples:", distinct_rank_psp_skin))
  print(paste("The proportion of identified ranks are:", (100*(1-proportion_na_psp_skin)), "with a total of NAs", na_count_psp_skin, "out of",total_taxa_psp_skin))
}

print("PSP URINE PELLET DATA:")
for(rank in taxonomic_ranks){
  rank_psp_urine<-taxa_psp_urine[, rank][!is.na(taxa_psp_urine)]
  distinct_rank_psp_urine<-length(unique(rank_psp_urine))
  total_taxa_psp_urine<-nrow(taxa_psp_urine)
  length(taxa_psp_urine[, rank][is.na(taxa_psp_urine[, rank])])
  na_count_psp_urine<-length(taxa_psp_urine[, rank][is.na(taxa_psp_urine[, rank])])
  proportion_na_psp_urine<-na_count_psp_urine/total_taxa_psp_urine
  print(paste("The number of distinct",rank,"(excluding NAs) in the PSP Urine samples:", distinct_rank_psp_urine))
  print(paste("The proportion of identified ranks are:", (100*(1-proportion_na_psp_urine)), "with a total of NAs", na_count_psp_urine, "out of",total_taxa_psp_urine))
}
print("BACT URINE PELLET DATA: ")
for(rank in taxonomic_ranks){
  rank_bact_urine<-taxa_bact_urine[, rank][!is.na(taxa_bact_urine)]
  distinct_rank_bact_urine<-length(unique(rank_bact_urine))
  total_taxa_bact_urine<-nrow(taxa_bact_urine)
  length(taxa_bact_urine[, rank][is.na(taxa_bact_urine[, rank])])
  na_count_bact_urine<-length(taxa_bact_urine[, rank][is.na(taxa_bact_urine[, rank])])
  proportion_na_bact_urine<-na_count_bact_urine/total_taxa_bact_urine
  print(paste("The number of distinct",rank,"(excluding NAs) in the Bact urine samples:", distinct_rank_bact_urine))
  print(paste("The proportion of identified ranks are:", (100*(1-proportion_na_bact_urine)), "with a total of NAs", na_count_bact_urine, "out of",total_taxa_bact_urine))
}

#NORMALIZATION
#CLR First
#Now I have three clean phyloseq objects per hypervariable region, each belongs to a different kit and thus extraction batch.
ps_psp_skin_clr<-norm_clr(psp_skin_decontaminated)
ps_psp_urine_clr<-norm_clr(psp_urine_decontaminated)
ps_bact_urine_clr<-norm_clr(bact_decontaminated)

#TSS next
#This method just performs normalization on the otu, so I need to recreate the ps objects
otu_psp_skin_tss<-norm_tss(otu_table(psp_skin_decontaminated))
otu_psp_urine_tss<-norm_tss(otu_table(psp_urine_decontaminated))
otu_bact_urine_tss<-norm_tss(otu_table(bact_decontaminated))

ps_psp_skin_tss<-phyloseq(otu_table(otu_psp_skin_tss, taxa_are_rows = FALSE), 
             sample_data(psp_skin_decontaminated),
             tax_table(psp_skin_decontaminated))
ps_psp_urine_tss<-phyloseq(otu_table(otu_psp_urine_tss, taxa_are_rows = FALSE),
              sample_data(psp_urine_decontaminated),
              tax_table(psp_urine_decontaminated))
ps_bact_urine_tss<-phyloseq(otu_table(otu_bact_urine_tss, taxa_are_rows = FALSE),
              sample_data(bact_decontaminated),
              tax_table(bact_decontaminated))

#TMM as last
otu_psp_skin_tmm<-norm_tmm(otu_table(psp_skin_decontaminated))
otu_psp_urine_tmm<-norm_tmm(otu_table(psp_urine_decontaminated))
otu_bact_urine_tmm<-norm_tmm(otu_table(bact_decontaminated))

ps_psp_skin_tmm<-phyloseq(otu_table(otu_psp_skin_tmm, taxa_are_rows = FALSE), 
                          sample_data(psp_skin_decontaminated),
                          tax_table(psp_skin_decontaminated))
ps_psp_urine_tmm<-phyloseq(otu_table(otu_psp_urine_tmm, taxa_are_rows = FALSE),
                           sample_data(psp_urine_decontaminated),
                           tax_table(psp_urine_decontaminated))
ps_bact_urine_tmm<-phyloseq(otu_table(otu_bact_urine_tmm, taxa_are_rows = FALSE),
                            sample_data(bact_decontaminated),
                            tax_table(bact_decontaminated))


#Calculating sum of ASVs per sample type
asv_total<-function(psp_skin_decontaminated, psp_urine_decontaminated, bact_decontaminated){
  df_psp_skin<-as.data.frame(sample_sums(psp_skin_decontaminated))
  df_psp_urine<-as.data.frame(sample_sums(psp_urine_decontaminated))
  df_bact_urine<-as.data.frame(sample_sums(bact_decontaminated))
  skin_metadata<-sample_data(psp_skin_decontaminated)
  urine_metadata<-sample_data(psp_urine_decontaminated)
  bact_metadata<-sample_data(bact_decontaminated)
  
  df_psp_skin$ID<-skin_metadata$ID
  df_psp_skin$sample_type<-skin_metadata$Type
  df_psp_skin$patient_id<-gsub(".*-(\\d{4})_.*", "\\1", df_psp_skin$ID)
  
  df_psp_urine$ID<-urine_metadata$ID
  df_psp_urine$sample_type<-urine_metadata$Type
  df_psp_urine$patient_id<-gsub(".*-(\\d{4})_.*", "\\1", df_psp_urine$ID)
  
  df_bact_urine$ID<-bact_metadata$ID
  df_bact_urine$sample_type<-bact_metadata$Type
  df_bact_urine$patient_id<-gsub(".*-(\\d{4})_.*", "\\1", df_bact_urine$ID)
  
  merged_dataframe_asv_count<-merge(df_psp_skin, df_psp_urine, by="patient_id",all=TRUE)
  merged_dataframe_asv_count<-merge(merged_dataframe_asv_count, df_bact_urine, by="patient_id", all=TRUE)
  #View(merged_dataframe_asv_count)
  write.csv(x=merged_dataframe_asv_count,"ASV_count_raw.csv")
}

asv_total(psp_skin_decontaminated, psp_urine_decontaminated, bact_decontaminated)

#Alpha diversity across non-normalized ps objects
alpha_estimation<-function(psp_skin_decontaminated, psp_urine_decontaminated, bact_decontaminated){
  alpha_diversity_psp_skin<-estimate_richness(psp_skin_decontaminated, measures=c("Shannon","Simpson","Observed","Chao1"))
  skin_metadata<-sample_data(psp_skin_decontaminated)
  alpha_diversity_psp_skin$sample_type<-skin_metadata$Type
  alpha_diversity_psp_skin$ID<-skin_metadata$ID
  alpha_diversity_psp_skin$patient_id<-gsub(".*-(\\d{4})_.*", "\\1", alpha_diversity_psp_skin$ID)

  alpha_diversity_psp_urine<-estimate_richness(psp_urine_decontaminated, measures=c("Shannon","Simpson","Observed","Chao1"))
  urine_metadata<-sample_data(psp_urine_decontaminated)
  alpha_diversity_psp_urine$sample_type<-urine_metadata$Type
  alpha_diversity_psp_urine$ID<-urine_metadata$ID
  alpha_diversity_psp_urine$patient_id<-gsub(".*-(\\d{4})_.*", "\\1", alpha_diversity_psp_urine$ID)
  alpha_diversity_psp_urine

  alpha_diversity_bact_urine<-estimate_richness(bact_decontaminated, measures=c("Shannon","Simpson","Observed","Chao1"))
  urine_bact_metadata<-sample_data(bact_decontaminated)
  alpha_diversity_bact_urine$sample_type<-urine_bact_metadata$Type
  alpha_diversity_bact_urine$ID<-urine_bact_metadata$ID
  alpha_diversity_bact_urine$patient_id<-gsub(".*-(\\d{4})_.*", "\\1", alpha_diversity_bact_urine$ID)
  alpha_diversity_bact_urine

  merged_dataframe<-merge(alpha_diversity_psp_skin, alpha_diversity_psp_urine, by="patient_id",all=TRUE)
  merged_dataframe<-merge(merged_dataframe, alpha_diversity_bact_urine, by="patient_id", all=TRUE)
  #View(merged_dataframe)
  write.csv(x=merged_dataframe,"richness_all_samples_non_normalized.csv")
}

alpha_estimation(psp_skin_decontaminated, psp_urine_decontaminated, bact_decontaminated)

#Rarefaction curves on non-normalized objects
ps_all_samples_raw<-merge_phyloseq(psp_skin_decontaminated, psp_urine_decontaminated, bact_decontaminated)
otu_all_samples_raw<-otu_table(ps_all_samples_raw)
otu_all_samples_raw
class(otu_all_samples_raw)<-"matrix" #Forcing it to accept it as matrix, as.matrix does not work
new_col_sample_kit<-data.frame(paste(sample_data(ps_all_samples_raw)$Kit, sample_data(ps_all_samples_raw)$Type2, sep="_"))
names(new_col_sample_kit)<-"Kit_and_sample"
new_col_sample_kit
sample_data(ps_all_samples_raw)$Kit_and_sample<-as.matrix(new_col_sample_kit)
rarefaction_data_raw_asv<-rarecurve(otu_all_samples_raw,step=100, label=FALSE
                ,col = as.factor(sample_data(ps_all_samples_raw)$Kit_and_sample
                ),tidy = TRUE)
rarefaction_data_raw_asv$patient_id<-gsub(".*-(\\d{4})_.*", "\\1", rarefaction_data_raw_asv$Site)
rarefaction_data_raw_asv$patient_id
rarefaction_df_asv_raw <- rarefaction_data_raw_asv %>%
pivot_wider(names_from = Sample, values_from = Species)
write.csv(rarefaction_df_asv_raw,"rarefaction_asv_count_raw.csv")

#BETA DIVERSITY ANALYSES
#Beta diversity analyses (cleaned and normalized)
ps_all_samples_clr<-merge_phyloseq(ps_psp_skin_clr, ps_psp_urine_clr, ps_bact_urine_clr)
ps_all_samples_tss<-merge_phyloseq(ps_psp_skin_tss, ps_psp_urine_tss, ps_bact_urine_tss)
ps_all_samples_tmm<-merge_phyloseq(ps_psp_skin_tmm, ps_psp_urine_tmm, ps_bact_urine_tmm)
#adding a new column that combines kit and sample type
ps_all_samples_clr
ps_all_samples_minus_empty_clr<-prune_samples(sample_sums(ps_all_samples_clr)>0, ps_all_samples_clr)
new_col_sample_kit<-data.frame(paste(sample_data(ps_all_samples_minus_empty_clr)$Kit, sample_data(ps_all_samples_minus_empty_clr)$Type2, sep="_"))
names(new_col_sample_kit)<-"Kit_and_sample"
sample_data(ps_all_samples_minus_empty_clr)$Kit_and_sample<-as.matrix(new_col_sample_kit)

ps_all_samples_tss
ps_all_samples_minus_empty_tss<-prune_samples(sample_sums(ps_all_samples_tss)>0, ps_all_samples_tss)
ps_all_samples_minus_empty_tss
sample_data(ps_all_samples_minus_empty_tss)
new_col_sample_kit<-data.frame(paste(sample_data(ps_all_samples_minus_empty_tss)$Kit, sample_data(ps_all_samples_minus_empty_tss)$Type2, sep="_"))
names(new_col_sample_kit)<-"Kit_and_sample"
sample_data(ps_all_samples_minus_empty_tss)$Kit_and_sample<-as.matrix(new_col_sample_kit)

ps_all_samples_tmm
ps_all_samples_minus_empty_tmm<-prune_samples(sample_sums(ps_all_samples_tmm)>0, ps_all_samples_tmm)
ps_all_samples_minus_empty_tmm
sample_data(ps_all_samples_minus_empty_tmm)
new_col_sample_kit<-data.frame(paste(sample_data(ps_all_samples_minus_empty_tmm)$Kit, sample_data(ps_all_samples_minus_empty_tmm)$Type2, sep="_"))
names(new_col_sample_kit)<-"Kit_and_sample"
sample_data(ps_all_samples_minus_empty_tmm)$Kit_and_sample<-as.matrix(new_col_sample_kit)


#Creating a function that gives me the beta plot and its permanova
graphing_beta_diversity<-function(ps_all_samples_minus_empty,method_plot,dist,string){
  getwd()
  sink(paste0(string,"_",method_plot,"_",dist,".txt"))
  group_colours<-c("deepskyblue","orange","blue")
  beta_dist_ord <- ordinate(ps_all_samples_minus_empty, method = method_plot, distance = dist)
  beta_dist_plot <- plot_ordination(ps_all_samples_minus_empty, ordination = beta_dist_ord, type = "samples", color = "Kit_and_sample")
  print(beta_dist_ord)
  plot <- beta_dist_plot + scale_color_manual(values=group_colours) + 
    stat_ellipse(alpha = 0.07, geom = "polygon", aes(fill = Kit_and_sample),show.legend=FALSE) +
    scale_fill_manual(values=group_colours) + 
    theme(
      panel.background = element_rect(fill = "white"),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "gray", size = 0.20),
      panel.grid.minor = element_line(colour = "gray", size = 0.05),
      panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
  distance_measure<-distance(ps_all_samples_minus_empty, method=dist)
  sampledf<-data.frame(sample_data(ps_all_samples_minus_empty))
  print(adonis2(distance_measure ~ Kit_and_sample, data=sampledf))
  print(adonis2(distance_measure ~ Kit, data=sampledf))
  print(adonis2(distance_measure ~ Type, data=sampledf))
  print(plot)
  sink()
}

#CLR Beta plots first
graphing_beta_diversity(ps_all_samples_minus_empty_clr,"NMDS","bray","clr")
ggsave("clr_NMDS_bray.pdf",width = 5.5,height = 3.5)
graphing_beta_diversity(ps_all_samples_minus_empty_clr,"NMDS","jaccard","clr")
ggsave("clr_NMDS_jaccard.pdf",width = 5.5,height = 3.5)
graphing_beta_diversity(ps_all_samples_minus_empty_clr,"PCoA","bray","clr")
ggsave("clr_pcoa_bray.pdf",width = 5.5,height = 3.5)
graphing_beta_diversity(ps_all_samples_minus_empty_clr,"PCoA","jaccard","clr")
ggsave("clr_pcoa_jaccard.pdf",width = 5.5,height = 3.5)

#TSS Beta plots next
graphing_beta_diversity(ps_all_samples_minus_empty_tss, "NMDS","bray","tss")
graphing_beta_diversity(ps_all_samples_minus_empty_tss, "NMDS","jaccard","tss")
graphing_beta_diversity(ps_all_samples_minus_empty_tss, "PCoA","bray","tss")
graphing_beta_diversity(ps_all_samples_minus_empty_tss, "PCoA","jaccard","tss")

#TMM Beta plots last
graphing_beta_diversity(ps_all_samples_minus_empty_tmm, "NMDS","bray","tmm")
graphing_beta_diversity(ps_all_samples_minus_empty_tmm, "NMDS","jaccard","tmm")
graphing_beta_diversity(ps_all_samples_minus_empty_tmm, "PCoA","bray","tss")
graphing_beta_diversity(ps_all_samples_minus_empty_tmm, "PCoA","jaccard","tmm")

#BAR PLOTS FOR THE GENERA AND SPECIFIC SAMPLE TYPES USING TRANSFORMED/NORMALIZED DATA
#This function can plot a relative abundance or absolute abundance plot
abundance_plots<-function(ps_object, top_to_look_for, type_of_plot, norm.method){
  clrs <- c("#e3869c","#5db746","#a25acd","#a3b334","#5d6bc8","#d29f36","#cf4ea6","#5bbe7d","#d4426c",
            "#4bb8aa","#cf4336","#5f96d0","#d96e2d","#bd8bd4","#487c3d","#974b7e","#a9af62","#a85147",
            "#7e6e28","#cc8957")
  file_name<-paste0(type_of_plot,"_",top_to_look_for,"_",norm.method)
  ps_object.genus<-tax_glom(ps_object, taxrank = "Genus")
  top<-names(sort(taxa_sums(ps_object.genus), decreasing=TRUE))[1:top_to_look_for]
  if(type_of_plot=="absolute"){
    ps_absolute_abundance<-prune_taxa(top, ps_object.genus)
    plot_genus_top<-plot_bar(ps_absolute_abundance, fill="Genus") + scale_fill_manual(values=clrs)
    plot_genus_top + geom_bar(aes(color=Genus, fill=Genus), stat="identity") + scale_color_manual(values=clrs) + theme(
      panel.background = element_rect(fill = "white"),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "gray", size = 0.20),
      panel.grid.minor = element_line(colour = "gray", size = 0.05)
    )
    write.csv(x=otu_table(ps_absolute_abundance), file = paste0(file_name,"_ASVcount.csv"))
    write.csv(x=tax_table(ps_absolute_abundance), file = paste0(file_name,"_taxaNames.csv"))
  }
  else if(type_of_plot=="relative"){
    # Calculate relative abundance at the genus level
    ps_relative_abundance <- transform_sample_counts(ps_object.genus, function(x) (x/sum(x)*100))
    ps_relative_abundance.genus<-prune_taxa(top, ps_relative_abundance)
    plot_genus_top <- plot_bar(ps_relative_abundance.genus, fill = "Genus") +
      scale_fill_manual(values = clrs) +
      geom_bar(aes(color = Genus, fill = Genus), stat = "identity") +
      scale_color_manual(values = clrs) +
      theme(
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "gray", size = 0.20),
        panel.grid.minor = element_line(colour = "gray", size = 0.05), box()
      )
    otu_df<-as.data.frame(otu_table(ps_relative_abundance.genus))
    otu_df<-otu_df %>%
      mutate(total_relative_abundance = rowSums(.))
    otu_df<- otu_df %>%
      mutate(new_column = 100 - total_relative_abundance)
    write.csv(x = otu_df, file=paste0(file_name,"_ASVcount.csv"))
    write.csv(x = tax_table(ps_relative_abundance.genus), file=paste0(file_name,"_taxaNames.csv"))
  }
  print(plot_genus_top)
  ggsave(filename = paste0(file_name,".pdf"),plot = plot_genus_top)
}

top_to_look_for=15

#absolute plots first
#clr
abundance_plots(ps_bact_urine_clr, top_to_look_for,"absolute","clr_bact_urine")
abundance_plots(ps_psp_skin_clr, top_to_look_for,"absolute","clr_psp_skin")
abundance_plots(ps_psp_urine_clr, top_to_look_for,"absolute","clr_psp_urine")
#tss
abundance_plots(ps_bact_urine_tss, top_to_look_for,"absolute","tss_bact_urine")
abundance_plots(ps_psp_skin_tss, top_to_look_for,"absolute","tss_psp_skin")
abundance_plots(ps_psp_urine_tss, top_to_look_for,"absolute","tss_psp_urine")
#tmm
abundance_plots(ps_bact_urine_tmm, top_to_look_for,"absolute","tmm_bact_urine")
abundance_plots(ps_psp_skin_tmm, top_to_look_for,"absolute","tmm_psp_skin")
abundance_plots(ps_psp_urine_tmm, top_to_look_for,"absolute","tmm_psp_urine")

#relative abundance plots next
abundance_plots(bact_decontaminated, top_to_look_for, "relative","bact_urine")
abundance_plots(psp_skin_decontaminated, top_to_look_for, "relative","psp_skin")
abundance_plots(psp_urine_decontaminated, top_to_look_for, "relative","psp_urine")

#Abundance plots using Maaslin2
ps_kit_object<-merge_phyloseq(bact_decontaminated,psp_urine_decontaminated)
ps_sample_object<-merge_phyloseq(psp_urine_decontaminated,psp_skin_decontaminated)
sample_data(ps_kit_object)

maaslin_function<-function(ps_object,type_to_analyze){
  metadata<-as.data.frame(sample_data(ps_all_samples_raw))
  class(metadata)<-"data.frame"
  ps_otu<-as(otu_table(ps_object), "matrix")
  ps_otu_df<-as.data.frame(ps_otu)
  fit_data = Maaslin2(input_data = ps_otu_df, 
                      input_metadata = metadata, 
                      output=type_to_analyze,
                      fixed_effects=c("Type2","Kit"))
  print(fit_data$results)
  #If age is available, perhaps this should be added as a random_effect into maaslin2
}
sample_data(ps_all_samples_raw)
maaslin_function(ps_all_samples_raw,"maaslin2_v1v2")
#View(tax_table(ps_all_samples_raw))

test2df#Saving a fasta file from the PS object first
ps_all_samples_raw %>%
  refseq() %>%
  Biostrings::writeXStringSet("ASV_raw.fasta")

fasta_file_unaligned<-read.fasta("ASV_raw.fasta")
entropy_values<-bio3d::entropy(fasta_file_unaligned)
entropy_values<-as.data.frame(entropy_values$H)
class(entropy_values)
window_size=10
entropy_values$rolling_average<-rollapply(entropy_values$`entropy_values$H`, width=window_size, FUN=mean, align="right", fill=NA)
entropy_values$rolling_sd<-rollapply(entropy_values$`entropy_values$H`, width=window_size, FUN=sd, align="right", fill=NA)
entropy_values
write.csv(entropy_values, "entropy_resolution.csv")

#Skin first
ps_psp_phyla_skin<-prune_samples(sample_data(ps_psp_phyla)$Type=="Penile-Sample", ps_psp_phyla)
ps_psp_phyla_skin_df<-as.data.frame(otu_table(ps_psp_phyla_skin))
ps_psp_phyla_skin_df_ordered<-rownames(ps_psp_phyla_skin_df)[order(rownames(ps_psp_phyla_skin_df))]
ps_psp_phyla_skin_df_ordered<-ps_psp_phyla_skin_df[ps_psp_phyla_skin_df_ordered,]
#Urine second
ps_psp_phyla_urine<-prune_samples(sample_data(ps_psp_phyla)$Type=="Urine-Sample", ps_psp_phyla)
ps_psp_phyla_urine_df<-as.data.frame(otu_table(ps_psp_phyla_urine))
ps_psp_phyla_urine_df_ordered<-rownames(ps_psp_phyla_urine_df)[order(rownames(ps_psp_phyla_urine_df))]
ps_psp_phyla_urine_df_ordered<-ps_psp_phyla_urine_df[ps_psp_phyla_urine_df_ordered,]
#Printing for prism graphing
cat(paste0(rownames(ps_psp_phyla_skin_df_ordered),"\t",ps_psp_phyla_skin_df_ordered$Firmicutes,"\n"))
cat(paste0(rownames(ps_psp_phyla_skin_df_ordered),"\t",ps_psp_phyla_skin_df_ordered$Proteobacteria,"\n"))
cat(paste0(rownames(ps_psp_phyla_skin_df_ordered),"\t",ps_psp_phyla_skin_df_ordered$Actinobacteriota,"\n"))
cat(paste0(rownames(ps_psp_phyla_skin_df_ordered),"\t",ps_psp_phyla_skin_df_ordered$Bacteroidota,"\n"))
cat(paste0(rownames(ps_psp_phyla_skin_df_ordered),"\t",ps_psp_phyla_skin_df_ordered$Campylobacterota,"\n"))
cat(paste0(rownames(ps_psp_phyla_skin_df_ordered),"\t",ps_psp_phyla_skin_df_ordered$Firmicutes_A,"\n"))
cat(paste0(rownames(ps_psp_phyla_skin_df_ordered),"\t",ps_psp_phyla_skin_df_ordered$Firmicutes_B,"\n"))
cat(paste0(rownames(ps_psp_phyla_skin_df_ordered),"\t",ps_psp_phyla_skin_df_ordered$Firmicutes_C,"\n"))
cat(paste0(rownames(ps_psp_phyla_skin_df_ordered),"\t",ps_psp_phyla_skin_df_ordered$Fusobacteriota,"\n"))
cat(paste0(rownames(ps_psp_phyla_skin_df_ordered),"\t",ps_psp_phyla_skin_df_ordered$Spirochaetota,"\n"))

cat(paste0(rownames(ps_psp_phyla_urine_df_ordered),"\t",ps_psp_phyla_urine_df_ordered$Firmicutes,"\n"))
cat(paste0(rownames(ps_psp_phyla_urine_df_ordered),"\t",ps_psp_phyla_urine_df_ordered$Proteobacteria,"\n"))
cat(paste0(rownames(ps_psp_phyla_urine_df_ordered),"\t",ps_psp_phyla_urine_df_ordered$Actinobacteriota,"\n"))
cat(paste0(rownames(ps_psp_phyla_urine_df_ordered),"\t",ps_psp_phyla_urine_df_ordered$Bacteroidota,"\n"))
cat(paste0(rownames(ps_psp_phyla_urine_df_ordered),"\t",ps_psp_phyla_urine_df_ordered$Campylobacterota,"\n"))
cat(paste0(rownames(ps_psp_phyla_urine_df_ordered),"\t",ps_psp_phyla_urine_df_ordered$Firmicutes_A,"\n"))
cat(paste0(rownames(ps_psp_phyla_urine_df_ordered),"\t",ps_psp_phyla_urine_df_ordered$Firmicutes_B,"\n"))
cat(paste0(rownames(ps_psp_phyla_urine_df_ordered),"\t",ps_psp_phyla_urine_df_ordered$Firmicutes_C,"\n"))
cat(paste0(rownames(ps_psp_phyla_urine_df_ordered),"\t",ps_psp_phyla_urine_df_ordered$Fusobacteriota,"\n"))
cat(paste0(rownames(ps_psp_phyla_urine_df_ordered),"\t",ps_psp_phyla_urine_df_ordered$Spirochaetota,"\n"))


ggplot2::p

