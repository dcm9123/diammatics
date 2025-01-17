#Daniel Castaneda Mogollon
#January 14th, 2025
#This script was created in order to run DADA2 against the w5w6 and w9w10 using
#my personally created 16S database using Sanger sequences

#Fethcing packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.18")
library("dada2")
library("devtools")


#Setting directory
path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/FemMicro16S/"
setwd(path)

#Importing non-chimeric reads from FemMicro16S
nonchimeric_w5 = readRDS("femmicro_output_w5w6/w5_seqtab_nochimeras.rds")
nonchimeric_w9 = readRDS("femmicro_output_w9w10/w8w9_seqtab_nochimeras.rds")
dim(nonchimeric_w5)
#Assigning taxonomy
assign_dada2_taxonomy<-function(nonchimeric_object, output_file){
  taxa <- assignTaxonomy(nonchimeric_object, refFasta = "../local_dada2_vsearch/dada2_new_updatad_taxonomy.fasta", 
                       multithread=TRUE,tryRC = TRUE,verbose = TRUE)
  asv_id<-paste0("ASV_", seq_along(colnames(nonchimeric_object)))
  taxonomy_df = as.data.frame(taxa)
  taxonomy_df$Sequence<-NULL
  taxonomy_df = cbind(ASV_ID = asv_id, taxonomy_df)
  head(taxonomy_df)
  taxa.print <- taxonomy_df # Removing sequence rownames for display only
  rownames(taxa.print) <- NULL
  #View(taxa.print)
  write.csv(x = taxa.print, file = output_file)
  #taxa.print
}

assign_dada2_taxonomy(nonchimeric_w5, "../local_dada2_vsearch/output/dada2_assignTaxonomy_w5w6_rdp.csv")
assign_dada2_taxonomy(nonchimeric_w9, "../local_dada2_vsearch/output/dada2_assignTaxonomy_w9w10_rdp.csv")
