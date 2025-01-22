#Daniel Castaneda Mogollon, PhD
#January 7th, 2025
#Script for Kat's project


#################################################################################
#INTRO###########################################################################
#################################################################################
library("phyloseq")
library("ggplot2")
library("pairwiseAdonis")
packageVersion("phyloseq")
path = "/Users/danielcm/Desktop/Sycuro/Projects/Chlamydia/"
setwd(path)
df = read.csv("vsearch_dada2_merged_clean_20241107_asv_collapse.csv", header=FALSE)
new_sample_names = as.character(unlist(df[1,66:455]))
asv_table = read.csv("kats_asv_table.csv")                                      #Reading the modified file for asv count and samples
asv_table_formatted = otu_table(asv_table, taxa_are_rows = TRUE)
asv_matrix = as.matrix(asv_table_formatted)                                     #Phyloseq needs a matrix to work properly, not a df
colnames(asv_matrix) = new_sample_names
asv_matrix
taxa_table = read.csv("kats_taxa_table.csv")                                    #Reading the modified file for taxonomy assignment per asv
taxa_matrix = as.matrix(taxa_table)
taxa_matrix = tax_table(taxa_matrix)
metadata = read.csv("kats_metadata_followup_cst.csv", row.names = 1)
metadata
#Establishing ps object
ps_object = phyloseq(asv_matrix,taxa_matrix)                                    #1173 taxa and 390 samples
ps_object = merge_phyloseq(ps_object, metadata)                                 #Sanity check
sample_data(ps_object)<-metadata
ps_object
#################################################################################
#DATA PREPROCESSING##############################################################
#################################################################################
ps_object_normalized = transform_sample_counts(ps_object, function(x) x / sum(x))
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

ps_filtered_1                               #66 taxon retained across 390 samples
ps_filtered_2                               #189 taxon retained across 390 samples
ps_filtered_3                               #413 taxon retained across 390 samples

ps1_temp = genefilter_sample(ps_filtered_1, filterfun_sample(function(x) x>0),   #This removes taxa that is not present in more than 10% of the samples
                             A=0.10*nsamples(ps_filtered_1))
ps2_temp = genefilter_sample(ps_filtered_1, filterfun_sample(function(x) x>0),   #This removes taxa that is not present in more than 5% of the samples
                             A=0.05*nsamples(ps_filtered_1))
ps1 = prune_taxa(ps1_temp, ps_filtered_1)                                        #30 taxa remain and 390 samples
ps2 = prune_taxa(ps2_temp, ps_filtered_1)                                        #44 taxa remain and 390 samples

ps1
ps2
#From now onwards, I will be using ps2, which has 0.001 cutoff mean of taxa across all samples, and the taxa is also present in 5% of the total samples


#################################################################################
#ALPHA DIVERSITY#################################################################
#################################################################################

#I will be calculating the alpha diversity on untrimmed data but removing non-existing ASVs per sample.
ps_non_zeros = prune_taxa(speciesSums(ps_object) > 0, ps_object)
alphas = estimate_richness(ps_non_zeros)
write.csv(x = alphas, file = "alpha_diversity_kat.csv")
alphas

#################################################################################
#BETA DIVERSITY##################################################################
#################################################################################

ps_css = transform_sample_counts(ps2, function(x) x/sum(x))
ps_css = transform_sample_counts(ps2, function(x) log1p(x))
otu_table(ps_css)                                                               #Sanity check for normalization 
sample_data(ps_css)
ps_css
sample_names(ps_css)

#Plot by urbanization
beta_plotting<-function(metadata_variable, dist, meth){
  ps_css_filtered = subset_samples(ps_css, !is.na(concat) & concat!="")
  beta_ordination = ordinate(ps_css_filtered, method = meth, distance = dist)
  group_colors = c("#5f9c9d","#d36f6f","#786a87")
  beta_plot = plot_ordination(ps_css, ordination = beta_ordination, type = "id2", color = metadata_variable)
  plot<-beta_plot + scale_color_manual(values = group_colors)+stat_ellipse(alpha = 0.07, geom = "polygon", aes(fill = metadata_variable), show.legend=FALSE) +
    theme(panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(colour = "gray", size=0.20),
          panel.grid.minor = element_line(colour = "gray", size = 0.05),
          panel.border = element_rect(colour="black", fill=NA, size = 1)
          )
  print(plot)
  filtered_sample_names = sample_names(ps_css_filtered)
  metadata_filtered = metadata[metadata$id2 %in% filtered_sample_names, ]
  distance_used = distance(ps_css_filtered, method = dist)
  print(pairwise.adonis2(distance_used ~ concat, data = metadata_filtered, method=dist,nperm = 999))
}

#Plot by
metadata
beta_plotting("concat","bray","PCoA")
beta_plotting("concat","jaccard","PCoA")
beta_plotting("VisitType","jaccard","PCoA")
beta_plotting("CST","bray","PCoA")
beta_plotting("CST","jaccard","PCoA")
beta_plotting("visitCST","bray","PCoA")
beta_plotting("visitCST","jaccard","PCoA")



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
