#Daniel Castaneda Mogollon, PhD
#January 9th, 2025
#This script was designed to use DESeq2 to compare the functional annotations
#found in the genomes of NS1, S2, S5, NS6 and determine if there is significant
#enrichment betweehn those and IAB+/-, as Anvi'o did not provide any significant
#findings (this is an alternative tool to do so).

library("DESeq2")
library("ggplot2")
path = ("/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/anvio_all_genomes_resequenced/02_CONTIGS/exported_tables/")
setwd(path)
df_merged_KOs = read.csv("merged_116_genomes_KOs.csv", header=TRUE)             #This file was merged using a Jupyter notebook script for KOs only
metadata_df = read.csv("genomes_metadata.csv", header=TRUE)
my_deseq2_function = function(df, metadata_file, output_title){
  deseq_object = DESeqDataSetFromMatrix(countData = df, 
                colData = metadata_file, design =~ consortia, tidy = TRUE)
  deseq_object
  deseq2_run <-DESeq(deseq_object)
  results_from_run = results(deseq2_run)
  print(results_from_run)
  write.csv(x = results_from_run,output_title)
}

my_deseq2_function(df_merged_KOs, metadata_df, "seroconversion_comparison.csv")
my_deseq2_function(df_merged_KOs, metadata_df, "consortia_comparison.csv")

metadata_df2 = read.csv("genomes_metadata_subsetns1s2.csv", header=TRUE)
df_merged_KOs_subset = read.csv("merged_116_genomes_KOs_subset_ns1s2.csv", header=TRUE)

my_deseq2_function(df_merged_KOs_subset, metadata_df2, "ns1_s2_deseq2_subset.csv")
