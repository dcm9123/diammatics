#Daniel Castaneda Mogollon, PhD
#January 7th, 2025
#Script for Kat's project


#################################################################################
#INTRO###########################################################################
#################################################################################
BiocManager::install("metagenomeSeq")
devtools::install_github("joey711/phyloseq")
library("zCompositions")
library("compositions")
library("phyloseq")
library("metagenomeSeq")
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
asv_matrix = as.matrix(asv_table_formatted)
asv_matrix#Phyloseq needs a matrix to work properly, not a df
colnames(asv_matrix) = new_sample_names
asv_matrix
taxa_table = read.csv("kats_taxa_table.csv")                                    #Reading the modified file for taxonomy assignment per asv
taxa_matrix = as.matrix(taxa_table)
taxa_matrix = tax_table(taxa_matrix)
metadata = read.csv("kats_metadata_cst_baseline.csv.csv", row.names = 1)
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
alphas = estimate_richness(ps_non_zeros)
write.csv(x = alphas, file = "alpha_diversity_kat.csv")
alphas

#################################################################################
#ALPHA DIVERSITY FROM SPECIES ###################################################
#################################################################################

df_alpha = read.csv("kats_alpha.csv", row.names=1)
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
  nsamples(ps_object)         #output is 116 for ps1
  min(sample_sums(ps_object))
  otu_mat <- as(otu_table(ps_object), "matrix")         # Convert OTU table to a matrix
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

otu_table(ps_css1) #sanity check
otu_table(ps_css2) #sanity check

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


#BETA-DIVERSITY PLOTTING
beta_plotting<-function(ps_object, metadata_variable, dist, meth, name){
  ps_filtered = subset_samples(ps_object, !is.na(concat) & concat!="")
  beta_ordination = ordinate(ps_filtered, method = meth, distance = dist)
  group_colors = c("#5f9c9d","#d36f6f","#786a87")
  beta_plot = plot_ordination(ps_object, ordination = beta_ordination, type = "id2", color = metadata_variable)
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
  metadata_filtered = metadata[metadata$id2 %in% filtered_sample_names, ]
  distance_used = distance(ps_filtered, method = dist)
  print(pairwise.adonis2(distance_used ~ concat, data = metadata_filtered, method=dist,nperm = 999))
}

#Plot by
metadata
beta_plotting(ps_css1, "concat","bray","PCoA", "pcoa_bray_cst_baseline_subset_001_and_10.tiff")
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
