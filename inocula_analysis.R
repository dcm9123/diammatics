#Daniel Castaneda Mogollon, PhD
#February 12th, 2025
#This script was generated to analyze the microbiome diversity of the inocula in T1D

library("phyloseq")
library("dada2")
library("metagenomeSeq")
library("zCompositions")
library("compositions")
library("ggplot2")
library("pairwiseAdonis")
library("Maaslin2")


path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/local_dada2_vsearch/"
setwd(path)

#################################################################################
#DATA PREPROCESSING##############################################################
#################################################################################

getting_inocula_asvs<-function(rds_file, taxa_file, plate_number){ #Need to fix the removal of SXX_L001 so the consortia won't be wrongly annotated
  df_inocula = readRDS(rds_file)
  sample_names = rownames(df_inocula)
  modified_sample_names = paste0(sample_names,"_plate",plate_number)
  rownames(df_inocula) = modified_sample_names
  names_of_files = sapply(strsplit(sample_names,"_L001"),`[`,1)
  consortia = logical(length(sample_names))
  seroconversion = logical(length(sample_names))
  inocula = logical(length(sample_names))
  id_and_plate = logical(length(sample_names))
  for (i in seq_along(sample_names)){
    if(grepl("_NS1_|NS1", sample_names[i])){
      consortia[i]="NS1"
      seroconversion[i]="IAB-"
    }
    else if(grepl("_S2_|S2", sample_names[i])){
      consortia[i]="S2"
      seroconversion[i]="IAB+"
    }
    else if(grepl("_NS6_|NS6", sample_names[i])){
      consortia[i]="NS6"
      seroconversion[i]="IAB-"
    }
    else if(grepl("_S5_|S5", sample_names[i])){
      consortia[i]="S5"
      seroconversion[i]="IAB+"
    }
    else{
      consortia[i]="Consortia wrong"
      seroconversion[i]="Seroconversion wrong"
    }
    if(grepl("_inocula_|_inoculum_|mouse", sample_names[i])){
      inocula[i] = 'Yes'
    }
    else{
      inocula[i] = 'No'
    }
    }
  #print(rds_file)
  metadata_generated = data.frame(ID=sample_names, Consortia=consortia, Seroconversion=seroconversion, Inocula=inocula)
  rownames(metadata_generated) = modified_sample_names
  taxa_table = read.table(taxa_file, header=TRUE, sep=',')
  #View(df_inocula)
  taxa_table = as.matrix(taxa_table)
  rownames(taxa_table)=taxa_table[,1] #This ensures that the sample names between taxa table and df_inocula match. This is very important! D:
  
  ps = phyloseq(otu_table(df_inocula, taxa_are_rows = FALSE), sample_data(metadata_generated), tax_table(taxa_table))
  dna = Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) = taxa_names(ps)
  ps = merge_phyloseq(ps, dna)
  return(ps)
}

ps_inocula1 = getting_inocula_asvs("plate1/seqtab_nochimeras.rds", "plate1/taxonomy/final_merged_tables/phyloseq_taxonomy.csv",1)
ps_inocula3 = getting_inocula_asvs("plate3/seqtab_nochimeras.rds", "plate3/taxonomy/final_merged_tables/phyloseq_taxonomy.csv",3)
ps_inocula4 = getting_inocula_asvs("plate4/seqtab_nochimeras.rds", "plate4/taxonomy/final_merged_tables/phyloseq_taxonomy.csv",4)
ps_inocula5 = getting_inocula_asvs("plate5/seqtab_nochimeras.rds", "plate5/taxonomy/final_merged_tables/phyloseq_taxonomy.csv",5)
sample_data(ps_inocula1)

ps_inocula1 #1346 ASVs and 38 samples and ctrls
ps_inocula3 #917 ASVs and 90 samples and ctrls
ps_inocula4 #1226 ASVs and 90 samples and ctls
ps_inocula5 #1442 ASVs and 56 samples and ctrls

merged_ps = merge_phyloseq(ps_inocula1, ps_inocula3, ps_inocula4, ps_inocula5) #The total of each ps alone adds up to 4,931 ASVs and 274 samples
sample_data(merged_ps)
#sample_data(merged_ps) Sanity check
merged_ps #This shows me 3,654 ASVs and 274 samples (previously 269), which suggests that there (before) are 7 repeated name IDs and 1,277 repeated ASV sequences during the merging
merged_ps_filtered = tax_glom(merged_ps, taxrank = "asv_seq", NArm = TRUE) #Sanity check, this makes sure that the number of ASVs do match the ones from my previous
                                                                           #merging process with 'merge_phyloseq'
tax_table(merged_ps) = tax_table(merged_ps)[,c(2:8,1)] #Had to do this as the first column is ASV seq and that is not what tax_table in ps expects, it expects kingdom first!
merged_ps_filtered = tax_glom(merged_ps, taxrank ='species_final')
merged_ps_filtered #This number went down to 81 taxa in total, suggesting the species merging worked!
sample_data(merged_ps_filtered)
merged_ps_inocula = subset_samples(physeq = merged_ps_filtered, Inocula!='No') #Saving only inocula
#sample_data(merged_ps_inocula) This gives me a total of 12 samples, which matches what I expect!
merged_ps_inocula_filtered = filter_taxa(merged_ps_inocula, function(x) sum(x)>0, prune=TRUE) #Removing taxa that do not have any counts
merged_ps_inocula_filtered #This number went down from 923 taxa to 518 or 60?

otu_table(merged_ps_inocula_filtered)
View(tax_table(merged_ps_inocula_filtered))

#################################################################################
#TAXA ANALYSIS###################################################################
#################################################################################
df_master = read.csv("/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/WGS_all_runs/FINAL_report_WGS021224.csv")
df_inocula_ns1 = subset(df_master, Selected.for.Downstream %in% "Yes" & Consortia %in% "NS1", select = c("Sample.ID", "Consortia", "GTDBtk.Species.Classification..pplacer."))
df_inocula_ns6 = subset(df_master, Selected.for.Downstream %in% "Yes" & Consortia %in% "NS6", select = c("Sample.ID", "Consortia", "GTDBtk.Species.Classification..pplacer."))
df_inocula_s2 = subset(df_master, Selected.for.Downstream %in% "Yes" & Consortia %in% "S2", select = c("Sample.ID", "Consortia", "GTDBtk.Species.Classification..pplacer."))
df_inocula_s5 = subset(df_master, Selected.for.Downstream %in% "Yes" & Consortia %in% "S5", select = c("Sample.ID", "Consortia", "GTDBtk.Species.Classification..pplacer."))

print(paste0("NS1 has a total of ",length(df_inocula_ns1$GTDBtk.Species.Classification..pplacer.), " strains, which add up to ",
             length(unique(df_inocula_ns1$GTDBtk.Species.Classification..pplacer.)), " species"))
print(paste0("NS6 has a total of ",length(df_inocula_ns6$GTDBtk.Species.Classification..pplacer.), " strains, which add up to ",
             length(unique(df_inocula_ns6$GTDBtk.Species.Classification..pplacer.)), " species"))
print(paste0("S2 has a total of ",length(df_inocula_s2$GTDBtk.Species.Classification..pplacer.), " strains, which add up to ",
             length(unique(df_inocula_s2$GTDBtk.Species.Classification..pplacer.)), " species"))
print(paste0("S5 has a total of ",length(df_inocula_s5$GTDBtk.Species.Classification..pplacer.), " strains, which add up to ",
             length(unique(df_inocula_s5$GTDBtk.Species.Classification..pplacer.)), " species"))


ns1_asv_genomes = unique(df_inocula_ns1$GTDBtk.Species.Classification..pplacer.)
ns6_asv_genomes = unique(df_inocula_ns6$GTDBtk.Species.Classification..pplacer.)
s2_asv_genomes = unique(df_inocula_s2$GTDBtk.Species.Classification..pplacer.)
s5_asv_genomes = unique(df_inocula_s5$GTDBtk.Species.Classification..pplacer.)

extracting_unique_species = function(ps_object){
  tax_table(ps_object)[,"species_final"] = paste(tax_table(ps_object)[, "genus_final"],
                                                 tax_table(ps_object)[, "species_final"],
                                                 sep=" ") #This merges genus and species into the species column
  tax_table(ps_object)[,"species_final"] = gsub("^([A-Za-z]+) \\1 ", "\\1 ", 
                                          tax_table(ps_object)[,"species_final"]) #This RE makes sure that if the first two strings are the same, it removes the 
  #first string (i.e. Akkermansia Akkermansia muciniphila -> Akkermansia muciniphila)
  tax_table(ps_object)[,"species_final"] <- gsub("^([A-Za-z]+ [A-Za-z]+).*", "\\1", tax_table(ps_object)[,"species_final"])
  ns1_asv_species_plate1 = subset_samples(physeq = ps_object, Consortia == "S5") #CHANGE THIS FOR EACH CONSORTIA!!!
  print(ns1_asv_species_plate1)
  ns1_asv_species_plate1 = filter_taxa(physeq = ns1_asv_species_plate1, flist = function(x) sum(x)>0, prune = TRUE)
  ns1_asv_species_plate1 = tax_glom(physeq = ns1_asv_species_plate1, taxrank = "species_final",NArm = TRUE)
  species_list = as.character(unique(tax_table(ns1_asv_species_plate1)[,"species_final"]))
  return(species_list)
}

#Getting the unique species from each plate and each consortia
ns1_plate1_species_list = extracting_unique_species(ps_inocula1) #Change for each inocula!
#none in ps_inocula3 ns1
ns1_plate4_species_list = extracting_unique_species(ps_inocula4)
ns1_plate5_species_list = extracting_unique_species(ps_inocula5)
ns6_plate1_species_list = extracting_unique_species(ps_inocula1)
ns6_plate3_species_list = extracting_unique_species(ps_inocula3)
ns6_plate4_species_list = extracting_unique_species(ps_inocula4)
#none in ps_inocula5 ns6
s2_plate1_species_list = extracting_unique_species(ps_inocula1)
s2_plate3_species_list = extracting_unique_species(ps_inocula3)
s2_plate4_species_list = extracting_unique_species(ps_inocula4)
s2_plate5_species_list = extracting_unique_species(ps_inocula5)

s5_plate1_species_list = extracting_unique_species(ps_inocula1)
s5_plate3_species_list = extracting_unique_species(ps_inocula3)
s5_plate4_species_list = extracting_unique_species(ps_inocula4)
#none in ps_inocula5 s5

overlap_and_differences = function(asv_genomes, plate_species){
  overlap_values = sort(intersect(asv_genomes, plate_species))
  print(paste0("The number of expected species is: ", length(asv_genomes)))
  expected_genus = sapply(strsplit(asv_genomes, " "),`[`,1)
  expected_genus = unique(sort(expected_genus))
  print(paste0("The number of expected genus are: ", length(expected_genus)))
  print(paste0("The genus expected are: ",paste(expected_genus, collapse=", ")))
  cat("\n")
  print(paste0("The overlap species are ", length(overlap_values),": ", paste(overlap_values, collapse=", ")))
  unique_to_genomes = sort(setdiff(asv_genomes, plate_species))
  print(paste0("The missed species to the genomes are ", length(unique_to_genomes), " :", paste(unique_to_genomes, collapse=", ")))
  unique_to_asvs = sort(setdiff(plate_species, asv_genomes))
  print(paste0("The unique species called by the asvs in the plates are: ", paste(unique_to_asvs, collapse=", ")))
}

overlap_and_differences(ns1_asv_genomes, ns1_plate1_species_list)
overlap_and_differences(ns1_asv_genomes, ns1_plate4_species_list)
overlap_and_differences(ns1_asv_genomes, ns1_plate5_species_list)
overlap_and_differences(ns6_asv_genomes, ns6_plate1_species_list)
overlap_and_differences(ns6_asv_genomes, ns6_plate3_species_list)
overlap_and_differences(ns6_asv_genomes, ns6_plate4_species_list)
overlap_and_differences(s2_asv_genomes, s2_plate1_species_list)
overlap_and_differences(s2_asv_genomes, s2_plate3_species_list)
overlap_and_differences(s2_asv_genomes, s2_plate4_species_list)
overlap_and_differences(s2_asv_genomes, s2_plate5_species_list)
overlap_and_differences(s5_asv_genomes, s5_plate1_species_list)
overlap_and_differences(s5_asv_genomes, s5_plate3_species_list)
overlap_and_differences(s5_asv_genomes, s5_plate4_species_list)


#################################################################################
#BETA DIVERSITY##################################################################
#################################################################################

#Normalizing by CSS
normalization_css<-function(ps_object){
  min(sample_sums(ps_object))
  #otu_mat <- as(otu_table(ps_object), "matrix")   # Convert OTU table to a matrix
  sample_nonzero_counts <- colSums(otu_table(ps_object) > 0)   # Count nonzero OTUs per sample
  summary(sample_nonzero_counts)                  # Check distribution of nonzero OTUs per sample
  ps_css1 = phyloseq_to_metagenomeSeq(ps_object)
  p1 = cumNormStat(ps_css1)
  #otu_table(ps1)
  ps_css1 = cumNorm(ps_css1, p = p1)
  ps_css1_ps = ps_object
  otu_table_ps_css1 = MRcounts(ps_css1,norm = TRUE)
  otu_table_ps_css1 = as.matrix(otu_table_ps_css1)
  print(otu_table_ps_css1)
#  otu_table_ps_css1
   taxa_are_rows_original = taxa_are_rows(ps_object)
   otu_table_ps_css1 = otu_table(otu_table_ps_css1, taxa_are_rows = taxa_are_rows_original)
#   print(otu_table_ps_css1)
#  otu_table(ps_css1_ps) = otu_table_ps_css1
   return(otu_table_ps_css1)
}

#merged_ps_inocula_filtered_m = phyloseq(
#  otu_table(otu_ps_inocula_css, taxa_are_rows=TRUE),
#  tax_table(merged_ps_inocula_filtered),
#  sample_data(merged_ps_inocula_filtered)
#)


normalization_clr<-function(ps_object){
  otu_matrix = as(otu_table(ps_object), "matrix")
  otu_matrix = cmultRepl(otu_matrix, label = 0, method = "CZM", z.warning = 0.99, z.delete = 0.99)
  otu_matrix_clr = t(apply(otu_matrix, 1, function(x){
    clr(x)
  }))
  otu_table_clr = otu_table(otu_matrix_clr, taxa_are_rows = taxa_are_rows(ps_object))
  ps_clr = ps_object
  otu_table(ps_clr)<-otu_table_clr
  return(ps_clr)
}

inocula_ps_clr = normalization_clr(merged_ps_inocula_filtered)

#BETA-DIVERSITY PLOTTING

beta_plotting<-function(ps_object, metadata_variable, dist, meth, name){
  beta_ordination = ordinate(ps_object, method = meth, distance = dist)
  #group_colors = c("#d36f6f","orange", "blue", "#5f9c9d") #change the number of colours if consortia vs seroconversion
  group_colors = c("#d36f6f", "#5f9c9d")
  beta_plot = plot_ordination(ps_object, ordination = beta_ordination, type = "samples", color = metadata_variable)
  plot<-beta_plot + scale_color_manual(values = group_colors)+
    stat_ellipse(alpha = 0.20, geom = "polygon", aes(fill = !!sym(metadata_variable)), show.legend=FALSE) +
    scale_fill_manual(values = group_colors) + 
    theme(panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(colour = "gray", size=0.20),
          panel.grid.minor = element_line(colour = "gray", size = 0.05),
          panel.border = element_rect(colour="black", fill=NA, size = 1)
    )
  print(plot)
  ggsave(filename = name,plot = plot,device = "pdf", units="in", width = 6, height=5, dpi=1200)
  metadata = as.data.frame(as.matrix(sample_data(ps_object))) #For some reason I need to save it as matrix, then as df for it to work
  rownames(metadata) = sample_names(ps_object)
  distance_used = distance(ps_object, method = dist)
  filtered_sample_names = sample_names(ps_object)
  #print(pairwise.adonis2(distance_used ~ Consortia, data = metadata, method=dist, nperm = 999))
  print(pairwise.adonis2(distance_used ~ Seroconversion, data = metadata, method=dist, nperm = 999)) #Change this for seroconversion or consortia
  
  #Change this for seroconversion or consortia
}


#Ploting the beta diversity
beta_plotting(inocula_ps_clr, "Consortia","euclidean","PCoA", "consortia.pdf")
beta_plotting(inocula_ps_clr, "Seroconversion", "euclidean","PCoA", "seroconversion.pdf")


#################################################################################
#ALPHA DIVERSITY & PAINTER PLOT##################################################
#################################################################################

#For alpha diversity we do not need to normalize, and we'll keep only the filter for
#no taxa with zeros

alphas = estimate_richness(merged_ps_inocula_filtered)
write.csv(x = alphas, file = "alpha_diversity_inocula.csv")
View(tax_table(merged_ps_inocula_filtered))


relative_abundance_plot = function(ps_object, top_to_look_for){
  print(rank_names(ps_object))
  ps_object_genus = tax_glom(ps_object, taxrank = "genus_final")
  print(ntaxa(ps_object_genus))
  top = names(sort(taxa_sums(ps_object_genus), decreasing=TRUE))[1:top_to_look_for]
  ps_relative_abundance = transform_sample_counts(ps_object_genus, function(x) (x/sum(x)*100))
  ps_relative_abundance_genus = prune_taxa(top, ps_relative_abundance)
  df = psmelt(ps_relative_abundance_genus)
  my_colors = c("#58643d","#6a42c3","#7ed05b","#c853bc","#cdb756",
               "#4b2e4b","#90c9b6","#c25d39","#8e92c2","#c35774")
  df$genus_final = factor(df$genus_final, levels = names(sort(tapply(df$Abundance, df$genus_final, sum), decreasing=TRUE)))
  unique_genera = levels(df$genus_final)
  repeated_colors = rep(my_colors, length.out = length(unique_genera))
  color_mapping = setNames(repeated_colors, unique_genera)
  
  plot_genus_top = ggplot(df, aes(x = Sample, y = Abundance, fill=genus_final)) +
    geom_bar(stat = "identity", position="stack") +
    theme_minimal() + labs(y = "Relative abundance (%)") + 
    theme(axis.text.x = element_text(angle=90)) +
    scale_fill_manual(values=color_mapping)
  print(plot_genus_top)
  #View(tax_table(ps_object_genus))
}
relative_abundance_plot(merged_ps_inocula_filtered, 31)



#Getting the .tsv file for picrust (not normalized)
file_tsv = t(otu_table(merged_ps_inocula_filtered))
file_tsv = as.data.frame(file_tsv)
file_tsv
write.table(file_tsv, file = "inocula_asv_file_12_samles.tsv", sep='\t')
write.csv(file_tsv, file = "inocula_asv_file_12_samples.tsv",sep='\t')

#Maaslin2 with the ASVs abundance generated table
###WORKING ON PICRUST2 DATA AND MAASLIN2
picrust_to_maaslin2<-function(path_picrust2, file_to_use, variable, normalization_method, transformation, ps, output_name){
  setwd(path_picrust2)
  #maaslin_pathways_df<-(read.delim("pathways_out/path_abun_unstrat_renamed.txt",header = TRUE))                            #I had to rename the file, as the pathways names were generating an error in Maaslin2
  maaslin_pathways_df<-read.delim(file_to_use)
  print(maaslin_pathways_df)
  #colnames(maaslin_pathways_df) <- sub("^.", "", colnames(maaslin_pathways_df), perl = TRUE)                                #R was adding an 'X' at the beginning of my sample names, so I removed that 'X'.
  maaslin_pathways_df<-t(maaslin_pathways_df)
  colnames(maaslin_pathways_df)<-maaslin_pathways_df[1,]
  maaslin_pathways_df<-maaslin_pathways_df[-1,]
  maaslin_pathways_df<-as.data.frame(maaslin_pathways_df)
  maaslin_pathways_df[] <- lapply(maaslin_pathways_df, as.numeric)
  metadata<-as.data.frame(sample_data(merged_ps_inocula_filtered))
  class(metadata)<-"data.frame"
  write.matrix(x=metadata,file = "metadata_file.txt")
  print(metadata)
  fit_data = Maaslin2(input_data = maaslin_pathways_df, input_metadata = metadata, output=output_name, 
                      fixed_effects=variable, normalization=normalization_method,transform = transformation)
}

picrust_to_maaslin2(path_picrust2 = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/local_dada2_vsearch/picrust2_inocula_merged_local_db/",
                    variable = "Seroconversion", normalization_method = "TSS", transformation = "LOG", file_to_use = "pred_metagenome_unstrat.tsv", ps=merged_ps_inocula_filtered,
                    output_name = 'maaslin2_results_localdb_tss_log2')

picrust_to_maaslin2(path_picrust2 = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/local_dada2_vsearch/picrust2_inocula_merged_default_db/KO_metagenome_out/",
                    variable = "Seroconversion", normalization_method = "TSS", transformation = "LOG", file_to_use = "pred_metagenome_unstrat.tsv", ps=merged_ps_inocula_filtered,
                    output_name = 'maaslin2_results_default_picrust_log2')
