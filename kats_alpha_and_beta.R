#Daniel Castaneda Mogollon, PhD
#January 7th, 2025
#Script for Kat's project
#R version = 4.4.3

#################################################################################
#INTRO###########################################################################
#################################################################################
R.version
#Installing
install.packages("BiocManager")
install.packages("remotes")
install.packages("devtools")
BiocManager::install("metagenomeSeq")
BiocManager::install("phyloseq", force = TRUE)
devtools::install_github("joey711/phyloseq")
install.packages("dendextend")
install.packages("ggdendro")
install.packages("zCompositions")
install.packages("compositions")
install.packages("pairwiseAdonis2")
install.packages("lattice")
BiocManager::install("pairwiseAdonis2")
BiocManager::install("biobakery/maaslin3")
devtools::install_github("joey711/phyloseq", force = TRUE, build_vignettes = FALSE)
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")


#Loading libraries
library("BiocManager")
library("zCompositions")
library("compositions")
library("phyloseq")
library("metagenomeSeq")
library("ggplot2")
library("pairwiseAdonis")
library("dendextend")
library("ggdendro")
library("dplyr")
library("maaslin3")

#Starting now
packageVersion("phyloseq")
path = "/Users/danielcm/Desktop/Sycuro/Projects/Chlamydia/"
setwd(path)
df = read.csv("vsearch_dada2_merged_clean_20241107_asv_collapse.csv", header=FALSE)
new_sample_names = as.character(unlist(df[1,66:455]))
asv_table = read.csv("kats_asv_table.csv")
asv_table_formatted = otu_table(asv_table, taxa_are_rows = TRUE)
asv_matrix = as.matrix(asv_table_formatted)                                                                     #Phyloseq needs a matrix to work properly, not a df
colnames(asv_matrix) = new_sample_names
taxa_table = read.csv("kats_taxa_table.csv")                                    #Here the ASVs are mapped to a specific taxa
taxa_matrix = as.matrix(taxa_table)
taxa_matrix = tax_table(taxa_matrix)
length(rownames(taxa_matrix))
metadata = read.csv("kats_metadata_confounders.csv", row.names = 1)

#Establishing ps object
ps_object = phyloseq(asv_matrix,taxa_matrix)                                    #1173 taxa and 390 samples
ps_object = merge_phyloseq(ps_object, metadata)                                 #Sanity check
ps_object = prune_samples(sample_names(ps_object)!="TKG20152702_S244_L001",ps_object)
otu_table(ps_object)
sample_data(ps_object) = metadata
ps_object
#View(df)

#################################################################################
#DATA PREPROCESSING##############################################################
#################################################################################
#Laura's way first#This calculates the relative abundance of each ASV sum across all samples, and discards if the threshold is
#below what I specify (0.001, 0.000001)

#To do: (i) Filter of 0.000001% and 5%, 0.000001% and nothing, and raw.
ps_non_zeros = prune_taxa(taxa_sums(ps_object) > 0, ps_object)                # Remove any ASV with 0 counts
ps_non_zeros
total_counts_per_asv = rowSums(otu_table(ps_non_zeros))                         # Add up all the ASV counts across all samples.
total_reads = sum(total_counts_per_asv)                                         # Get the total number of reads across all ASVs and all samples
relative_abundance = total_counts_per_asv/total_reads                           # Calculate the relative abundance for each ASV
relative_abundance                                                              # Sanity check
thresholds = c(0.001, 0.00001)
ps_asv_filter= prune_taxa(relative_abundance>=thresholds[2], ps_non_zeros)
#laura_ps2 = prune_taxa(relative_abundance>=thresholds[2], ps_non_zeros)                                                                  #Down to 185 taxa
ps_temp = genefilter_sample(ps_asv_filter, filterfun_sample(function(x) x>0),
                                             A = 0.05*nsamples(ps_asv_filter))
ps_asv_and_sample_filter = prune_taxa(ps_temp, ps_asv_filter) 
ps3 = ps_asv_and_sample_filter #ps1 to use (hard filter)
ps2 = ps_asv_filter #ps2 to use (relaxed filter)
ps1 = ps_non_zeros #ps3 to use (raw)

ps1 #501 ASVs
ps2 #400 ASVs
ps3 #80 ASVs
ps_object #1173 taxa (original with no filtering) #Skip line 92 to line 183 if not doing Alpha-diversity


alphas1 = phyloseq::estimate_richness(ps1) #Not typically used
alphas2 = phyloseq::estimate_richness(ps2) #Not typically used
alphas3 = phyloseq::estimate_richness(ps3) #Not typically used
alphas4 = phyloseq::estimate_richness(ps_object) #Use unfiltered one!

alphas4 = cbind(id2 = rownames(alphas4), alphas4)
alphas4 = merge(alphas4, metadata[,c("id2","VisitType")], by="id2", all.x=TRUE)
alphas4 #Sanity check

alphas1 = cbind(id2 = rownames(alphas1), alphas1)
alphas1 = cbind(id2 = rownames(alphas1), alphas1)
alphas1
metadata
alphas1 = merge(alphas1, metadata[,c("id2","VisitType")], by="id2", all.x=TRUE)
write.csv(alphas1, "alphas_no_zeroes.csv")

otu_values = otu_table(ps_object, taxa_are_rows = TRUE)
otu_table(ps_object)
total_abundance = colSums(otu_values)
max_abundance = apply(otu_values,2,max)
berger_index = max_abundance/total_abundance
sample_data(ps_object)$berger = as.character(berger_index)
sample_data(ps_object)
x <- as.matrix(sample_data(ps_object))

write.csv(x, "alpha_raw_with_berger.csv")
x = as.data.frame(x)
class(x)
setdiff(sample_data(ps_object)$id2, as.character(alphas4$id2))
dim(x)

################################################################################
################################################################################
################################################################################
#Daniel's way second
ps_non_zeros = prune_taxa(speciesSums(ps_object) > 0, ps_object)
ps_object_normalized = transform_sample_counts(ps_non_zeros, function(x) x / sum(x))
thresholds = c(0.001,0.0001,0.00001)
filtered_ps = lapply(thresholds, function(thres){
  keep_taxa = taxa_names(ps_object_normalized)[
    apply(otu_table(ps_object_normalized), 1, function (x) mean(x) > thres)
  ]
  prune_taxa(keep_taxa, ps_object)
})
ps_filtered_1 = filtered_ps[[1]]
ps_filtered_2 = filtered_ps[[2]]
ps_filtered_3 = filtered_ps[[3]]

ps_filtered_1                               #66 taxon retained across 117 samples
ps_filtered_2                               #189 taxon retained across 117 samples
ps_filtered_3                               #413 taxon retained across 117 samples

ps1_temp = genefilter_sample(ps_filtered_1, filterfun_sample(function(x) x>0),   #This removes taxa that is not present in more than 10% of the samples
                             A=0.10*nsamples(ps_filtered_1))
ps2_temp = genefilter_sample(ps_filtered_3, filterfun_sample(function(x) x>0),   #This removes taxa that is not present in more than 5% of the samples
                             A=0.05*nsamples(ps_filtered_1))
ps1 = prune_taxa(ps1_temp, ps_filtered_1)                                        #30 taxa remain and 390 samples 0.001 and 10%
ps2 = prune_taxa(ps2_temp, ps_filtered_3)                                        #44 taxa remain and 390 samples 0.00001 and 5%

ps1 #0.001 and 10%
ps2 #0.00001 and 5%


#################################################################################
#ALPHA DIVERSITY FROM ASVS#######################################################
#################################################################################

#I will be calculating the alpha diversity on untrimmed data but removing non-existing ASVs per sample.
ps_non_zeros = prune_taxa(speciesSums(ps_object) > 0, ps_object)
ps_non_zeros
alphas = estimate_richness(ps_non_zeros)
write.csv(x = alphas, file = "alpha_diversity_kat.csv")
alphas

#################################################################################
#ALPHA DIVERSITY FROM SPECIES ###################################################
#################################################################################

df_alpha = read.csv("kats_alpha.csv", row.names=1)
df_alpha
View(df_alpha)
df_alpha = as.matrix(df_alpha)
df_alpha = apply(df_alpha, 2, as.numeric)
df_alpha
species_table = otu_table(df_alpha, taxa_are_rows = TRUE)
ps_alpha = phyloseq(species_table)
ps_alpha_non_zeroes = prune_taxa(taxa_sums(ps_alpha)>0, ps_alpha)
alphas_species = estimate_richness(ps_alpha)
alphas_species = write.csv(x = alphas_species,file = "alpha_diversity_kat.csv")

#################################################################################
#BETA DIVERSITY##################################################################
#################################################################################

#Normalizing by CSS
normalization_css<-function(ps_object){
  ntaxa(ps_object)            #output is 30 for ps1
  nsamples(ps_object)         #output is 117 for ps1
  min(sample_sums(ps_object))
  otu_mat <- as(otu_table(ps_object), "matrix")   # Convert OTU table to a matrix
  sample_nonzero_counts <- colSums(otu_mat > 0)   # Count nonzero OTUs per sample
  summary(sample_nonzero_counts)                  # Check distribution of nonzero OTUs per sample
  ps_css1 = phyloseq_to_metagenomeSeq(ps_object)
  p1 = cumNormStat(ps_css1)
  #otu_table(ps1)
  ps_css1 = cumNorm(ps_css1, p = p1)
  ps_css1_ps = ps_object
  otu_table_ps_css1 = MRcounts(ps_css1,norm = TRUE)
  otu_table_ps_css1 = as.matrix(otu_table_ps_css1)
  otu_table_ps_css1
  taxa_are_rows_original = taxa_are_rows(ps_object)
  otu_table_ps_css1 = otu_table(otu_table_ps_css1, taxa_are_rows = taxa_are_rows_original)
  otu_table(ps_css1_ps) = otu_table_ps_css1
  return(ps_css1_ps)
}

ps_css1 = normalization_css(ps1)
ps_css2 = normalization_css(ps2)
ps_css3 = normalization_css(ps3)

otu_table(ps_css1) #sanity check
otu_table(ps_css2) #sanity check
otu_table(ps_css3) #sanity check

#Normalizing by CLR
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
  #otu_table(ps_clr)
}

ps_clr1 = normalization_clr(ps1)
ps_clr2 = normalization_clr(ps2)
ps_clr1
ps_clr2

otu_table(ps_clr2)

#EUCLIDEAN CLUSTERING
euclidean_dendrogram = function(ps_object){
  euclidean_dist = dist(otu_table(t(ps_object)), method="euclidean")
  clustering = hclust(euclidean_dist, method="complete")
  cst_type = sample_data(ps_object)[["CST"]]
  names(cst_type) = sample_names(ps_object)
  unique_values = unique(cst_type[cst_type!=""])
  colors = setNames(rainbow(length(unique_values)), unique_values)
  colors[""] = "black"
  colors["CSTI"]="blue"
  colors["CSTIII"]="purple"
  colors["CSTIV"]="firebrick1"
  ordered_samples = rownames(t(otu_table(ps_object))[clustering$order])
  for(sample in ordered_samples){
    cat(sample,sep = "\n")
  }
  label_colors = colors[cst_type[ordered_samples]]
  
  # Adjust margins to create more space for the text below the plot
  par(mfrow = c(1,1), xpd=NA)  # Increase the bottom margin for text labels
  
  # Plot dendrogram without labels, but shrink the height to allow space below
  plot(clustering, labels = FALSE, cex = 0.5, hang=0.2,
       main = "", xlab = "", sub = "", 
       ylim = c(-max(clustering$height) * 10, max(clustering$height) * 1.2))  # Adjust ylim
  
  # Get leaf x-positions for text labels
  dend = as.dendrogram(clustering)
  leaf_x_positions = order.dendrogram(dend)
  # Add text labels below the dendrogram
  # Adjust the 'y' position to move text further down
  text(leaf_x_positions, 
       labels = ordered_samples, 
       col = label_colors, 
       srt = 90, adj = 1.5, cex = 0.5)
}
ps_css1_baseline = subset_samples(ps_css1,VisitType=="Baseline")
ps_css2_baseline = subset_samples(ps_css2,VisitType=="Baseline")
ps_css3_baseline = subset_samples(ps_css3,VisitType=="Baseline")
ps_css1_followup = subset_samples(ps_css1,VisitType=="Follow Up")
ps_css2_followup = subset_samples(ps_css2,VisitType=="Follow Up")
ps_css3_followup = subset_samples(ps_css3,VisitType=="Follow Up")
ps1_baseline = subset_samples(ps1,VisitType=="Baseline")
ps2_baseline = subset_samples(ps2,VisitType=="Baseline")
ps3_baseline = subset_samples(ps3,VisitType=="Baseline")
ps1_followup = subset_samples(ps1,VisitType=="Follow Up")
ps2_followup = subset_samples(ps2,VisitType=="Follow Up")
ps3_followup = subset_samples(ps3,VisitType=="Follow Up")

euclidean_dendrogram(ps_css3_baseline)
euclidean_dendrogram(ps_css2_baseline)
euclidean_dendrogram(ps_css1_baseline)
euclidean_dendrogram(ps_css3_followup)
euclidean_dendrogram(ps_css2_followup)
euclidean_dendrogram(ps_css1_followup)

euclidean_dendrogram(ps3_baseline)
euclidean_dendrogram(ps2_baseline)
euclidean_dendrogram(ps1_baseline)
euclidean_dendrogram(ps3_followup)
euclidean_dendrogram(ps2_followup)
euclidean_dendrogram(ps1_followup)


#BETA-DIVERSITY PLOTTING
beta_plotting<-function(ps_object, metadata_variable, dist, meth, name){
  ps_filtered = subset_samples(ps_object, !is.na(Yrs_quartile) & Yrs_quartile!="")
  beta_ordination = ordinate(ps_filtered, method = meth, distance = "bray")
  group_colors = c("#5f9c9d","#d36f6f","#786a87","#c4b7a6")
  beta_plot = plot_ordination(ps_filtered, ordination = beta_ordination, type = "samples", color = metadata_variable)
  plot<-beta_plot + scale_color_manual(values = group_colors)+stat_ellipse(alpha = 0.20, geom = "polygon", aes(fill = !!sym(metadata_variable)), show.legend=FALSE) +
    scale_fill_manual(values = group_colors) + 
    theme(panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(colour = "gray", size=0.20),
          panel.grid.minor = element_line(colour = "gray", size = 0.05),
          panel.border = element_rect(colour="black", fill=NA, size = 1)
          )
  print(plot)
  ggsave(filename = name,plot = plot,device = "tiff", units="in", width = 6, height=5, dpi=1200)
  filtered_sample_names = sample_names(ps_filtered)
  metadata_filtered = metadata[metadata$id2 %in% filtered_sample_names,]
  distance_used = phyloseq::distance(physeq = ps_filtered,"bray")
  adonis_results = pairwise.adonis2(distance_used ~ Yrs_quartile, data = metadata_filtered, method="bray",nperm = 999,)
  adonis_results$p.adjusted = p.adjust(adonis_results$p.value, method = "fdr")
  print(adonis_results)
}

#Plot by


beta_plotting(ps_css1, "concat","bray","NMDS", "NMDS_bray_cst4.tiff")
beta_plotting(ps_css2, "concat","bray","NMDS", "NMDS_bray_cst4_soft.tiff")
beta_plotting(ps_css3, "concat","bray","NMDS", "NMDS_bray_cst4_hard.tiff")
beta_plotting(ps_css1, "concat","jaccard","NMDS", "NMDS_jaccard_cst4.tiff")
beta_plotting(ps_css2, "concat","jaccard","NMDS", "NMDS_jaccard_cst4_soft.tiff")
beta_plotting(ps_css3, "concat","jaccard","NMDS", "NMDS_jaccard_cst4_hard.tiff")

beta_plotting(ps_css1, "concat","jaccard","PCoA","pcoa_jaccard_cst_baseline_subset_001_and_10.tiff")
beta_plotting(ps_css2, "concat","bray","PCoA","pcoa_bray_cst_baseline_subset_00001_and_5.tiff")
beta_plotting(ps_css2, "concat","jaccard","PCoA","pcoa_jaccard_cst_baseline_subset_00001_and_5.tiff")
beta_plotting(ps_css1, "concat", "euclidean","PCoA","PCA_euclidean_cst_baseline_subset_001_and_10.tiff")
beta_plotting(ps_css2, "concat", "euclidean", "PCoA", "PCA_eulcidean_cst_baseline_subset_00001_and_5.tiff")
beta_plotting(ps_clr1, "concat", "euclidean", "PCoA","pcoa_euclidean_cst_baseline_subset_clr_001_and_10.tiff")
beta_plotting(ps_clr2, "concat", "euclidean", "PCoA", "pcoa_euclidean_cst_baseline_subset_clr_00001_and_5.tiff")

beta_plotting(ps_css1, "concat","bray","PCoA", "pcoa_bray_cst_follow_subset_001_and_10.tiff")
beta_plotting(ps_css1, "concat","jaccard","PCoA","pcoa_jaccard_cst_follow_subset_001_and_10.tiff")
beta_plotting(ps_css2, "concat","bray","PCoA","pcoa_bray_cst_follow_follow_00001_and_5.tiff")
beta_plotting(ps_css2, "concat","jaccard","PCoA","pcoa_jaccard_cst_follow_follow_00001_and_5.tiff")
beta_plotting(ps_css1, "concat", "euclidean","PCoA","PCA_euclidean_cst_follow_001_and_10.tiff")
beta_plotting(ps_css2, "concat", "euclidean", "PCoA", "PCA_eulcidean_cst_follow_subset_00001_and_5.tiff")
beta_plotting(ps_clr1, "concat", "euclidean", "PCoA","pcoa_euclidean_cst_follow_subset_clr_001_and_10.tiff")
beta_plotting(ps_clr2, "concat", "euclidean", "PCoA", "pcoa_euclidean_cst_follow_subset_clr_00001_and_5.tiff")

beta_plotting(ps_css1, "had_debut","bray","NMDS", "NMDS_bray_debut.tiff")
beta_plotting(ps_css2, "had_debut","bray","NMDS", "NMDS_bray_debut_soft.tiff")
beta_plotting(ps_css3, "had_debut","bray","NMDS", "NMDS_bray_debut_hard.tiff")

beta_plotting(ps_css1, "condom_use","bray","NMDS", "NMDS_bray_condom.tiff")
beta_plotting(ps_css2, "condom_use","bray","NMDS", "NMDS_bray_condom_soft.tiff")
beta_plotting(ps_css3, "condom_use","bray","NMDS", "NMDS_bray_condom_hard.tiff")

beta_plotting(ps_css1, "Yrs_quartile","bray","NMDS", "NMDS_bray_yrs_debut.tiff")
beta_plotting(ps_css2, "Yrs_quartile","bray","NMDS", "NMDS_bray_yrs_debut_soft.tiff")
beta_plotting(ps_css3, "Yrs_quartile","bray","NMDS", "NMDS_bray_yrs_debut_hard.tiff")



metadata$Yrs_since_debut

sample_data(ps_css1)

pairwiseAdonis::


ps_clr1
ps_clr2

beta_plotting("urban","jaccard","PCoA")
beta_plotting("housing_condition","bray","NMDS")
beta_plotting("age_at_debut","bray","NMDS")
beta_plotting("age_at_menarche","bray","NMDS")
beta_plotting("monthly_income","bray","NMDS")
beta_plotting("hormone","bray","NMDS")
beta_plotting("urban","jaccard","NMDS")
beta_plotting("housing_condition","jaccard","NMDS")
beta_plotting("age_at_debut","jaccard","NMDS")
beta_plotting("age_at_menarche","jaccard","NMDS")
beta_plotting("monthly_income","jaccard","NMDS")
beta_plotting("hormone","jaccard","NMDS")

#################################################################################
#DIFFERENTIAL ABUNDANCE ANALYSIS#################################################
#################################################################################
ps_object #Use the raw object without any filtering or removal
sum(is.na(as.data.frame(tax_table(ps_object))$Species)) #358 are NAs in the Species column
sum(is.na(as.data.frame(tax_table(ps_object))$Genus)) #Only 41 are NAs in the Genus column

agglomeration = function(ps_object, with_NAs){
  original_otu_table = as.data.frame(otu_table(ps_object)) #Retains the original otu_table
  original_taxa_table = as.data.frame(tax_table(ps_object))#Retains the original taxa table
  removed_x_taxa_table = original_taxa_table[,colnames(original_taxa_table)!="X"]
  removed_x_taxa_table$genus_species = paste(removed_x_taxa_table$Genus,removed_x_taxa_table$Species)
  removed_x_taxa_table$genus_species = gsub("\\s+"," ",removed_x_taxa_table$genus_species) #Ensures I only have one white space between genus and species
  removed_x_taxa_table = removed_x_taxa_table[removed_x_taxa_table$genus_species!="NA NA",]
  tax_table(ps_object) = as.matrix(removed_x_taxa_table)
  if(with_NAs==TRUE){
    ps_object = tax_glom(ps_object,taxrank = "genus_species", NArm=FALSE)
  }
  else{
    ps_object = tax_glom(ps_object,taxrank = "genus_species", NArm = TRUE)
  }
  df_tax_table = as.data.frame(tax_table(ps_object))
  row_names_taxa_table = df_tax_table$genus_species
  tax_table(ps_object) = NULL
  row.names(otu_table(ps_object)) = row_names_taxa_table
  row.names(df_tax_table) = row_names_taxa_table
  tax_table(ps_object) = as.matrix(df_tax_table)
  unique_taxa = length(unique(tax_table(ps_object)[,"genus_species"])) #Counts the number of unique genus species denomination
  length_taxa = length(tax_table(ps_object)[,"Species"]) #Counts the number of rows in the species column
  print(paste0("The number of unique species is ",unique_taxa,
               " and the number of total species is ",length_taxa)) #This numbers should match!
  return(ps_object)
}

ps_species = agglomeration(ps_object,with_NAs = FALSE)
ps_species #From ~1200 to 510 taxa after removing NAs and grouping same species!
sample_data(ps_species)

#Pruning bps_object#Pruning by visit type empty
visit_type = as.character(sample_data(ps_species)$VisitType)
good_visit = !is.na(visit_type) & !grepl("^[[:space:]]*$", visit_type) #Removes empty spaces into account 
ps_species_visit = prune_samples(good_visit, ps_species)
sample_data(ps_species_visit)

#Pruning by CST empty
cst_type = as.character(sample_data(ps_species)$CST)
good_cst = !is.na(cst_type) & !grepl("^[[:space:]]*$", cst_type) #Removes empty spaces into account 
ps_species_cst = prune_samples(good_cst, ps_species)

#Pruning by Chlamydia VisitType & zeros
ps_species_lacto = subset_samples(physeq = ps_species,Lacto_Dominance=="Dom") #From 116 samples to 83
ps_species_lacto = subset_samples(ps_species_lacto,VisitType!="Chlamydia")#Went down to 54 samples
ps_species_lacto = prune_taxa(taxa_sums(ps_species_lacto) > 0, ps_species_lacto) #Species from 510 to 116

ps_species_cst1 = subset_samples(ps_species,CST=="CSTI") #From 116 samples to 33
ps_species_cst1 = subset_samples(ps_species_cst1, VisitType!="Chlamydia") #Went down to 24
ps_species_cst1 = prune_taxa(taxa_sums(ps_species_cst1)>0, ps_species_cst1) #Species from 510 to 107 taxa


#Running maaslin3, only species abundance and metadata file are required
metadata_species_visit = as.data.frame(as.matrix(sample_data(ps_species_visit)))
otu_table_species_visit = t(as.data.frame(otu_table(ps_species_visit)))
metadata_species_cst = as.data.frame(as.matrix(sample_data(ps_species_cst)))
otu_table_species_cst = t(as.data.frame(otu_table(ps_species_cst)))
metadata_species_lacto = as.data.frame(as.matrix(sample_data(ps_species_lacto)))
otu_table_species_lacto = t(as.data.frame(otu_table(ps_species_lacto)))
metadata_species_cst1 = as.data.frame(as.matrix(sample_data(ps_species_cst1)))
otu_table_species_cst1 = t(as.data.frame(otu_table(ps_species_cst1)))


maaslin3(input_data = otu_table_species_visit,input_metadata = metadata_species_visit,
         output = "Visit_type_maaslin3",fixed_effects = "VisitType",
         normalization = "TSS",transform = "LOG")

maaslin3(input_data = otu_table_species_cst, input_metadata = metadata_species_cst, 
         output = "CST1_maaslin3", fixed_effects = "CST", normalization = "TSS",
         transform = "LOG", reference = "CST,CSTI")

maaslin3(input_data = otu_table_species_cst, input_metadata = metadata_species_cst, 
         output = "CST3_maaslin3", fixed_effects = "CST", normalization = "TSS",
         transform = "LOG", reference = "CST,CSTIII")

maaslin3(input_data = otu_table_species_lacto, input_metadata = metadata_species_lacto,
         output = "Lacto_dominant_raw_maaslin3", fixed_effects = "VisitType", normalization = "TSS",
         transform = "LOG")

maaslin3(input_data = otu_table_species_lacto, input_metadata = metadata_species_lacto,
         output = "Lacto_dominant_soft_maaslin3", fixed_effects = "VisitType", normalization = "TSS",
         transform = "LOG",min_abundance = 0.000001) #76 species after filtering

maaslin3(input_data = otu_table_species_lacto, input_metadata = metadata_species_lacto,
         output = "Lacto_dominant_hard_maaslin3", fixed_effects = "VisitType", normalization = "TSS",
         transform = "LOG",min_abundance = 0.000001, min_prevalence = 0.05) 

maaslin3(input_data = otu_table_species_cst1, input_metadata = metadata_species_cst1,
         output = "CSTI_raw_maaslin3", fixed_effects = "VisitType", normalization = "TSS",
         transform = "LOG")

maaslin3(input_data = otu_table_species_cst1, input_metadata = metadata_species_cst1,
         output = "CSTI_soft_maaslin3", fixed_effects = "VisitType", normalization = "TSS",
         transform = "LOG", min_abundance = 0.000001)

maaslin3(input_data = otu_table_species_cst1, input_metadata = metadata_species_cst1,
         output = "CSTI_hard_maaslin3", fixed_effects = "VisitType", normalization = "TSS",
         transform = "LOG", min_abundance = 0.000001, min_prevalence = 0.05)

sample_data(ps_species)
