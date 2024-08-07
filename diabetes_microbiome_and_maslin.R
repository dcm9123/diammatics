###FEMMICRO PHYLOSEQ ANALYSES###
#Daniel Castaneda Mogollon
#July 9th, 2024

#Installing packages
install.packages('devtools')
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

remove.packages("Maaslin2")
purge("Maaslin2")
remotes::install_github(repo="biobakery/Maaslin2", force=TRUE)
remotes::install_github("biobakery/Maaslin2")

#Loading packages
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



###INITIALIZING PATH AND FILE LOCATION
path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/FemMicro16S/femmicro_output_w8w9"
setwd(path)                                                                     #Setting the path of where Im working

###INITIALIZING METADATA FILE STRUCTURING
seqtab.nochim<-readRDS(file = "w8w9_seqtab_nochimeras.rds")                     #Reading the non-chimeric sequences from FemMicro
samples.out<-rownames(seqtab.nochim)                                            #Getting the name of the samples
samples.out
names_of_files <- sapply(strsplit(samples.out,"_L001"),`[`,1)                   #This divides the sample names by week, sex, and id
mice_sex<-logical(length(names_of_files))                                       #Creating a variable to store the sex of the mice
community_type<-logical(length(names_of_files))                                 #Creating a variable to store the community type (IAB+/-)
community_subtype<-logical(length(names_of_files))                              #Creating a variable to store the community subtype (NS1,NS6,S2,S5)
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
samdf<-data.frame(ID=samples.out, Sex=mice_sex, Community=community_type, 
                  Subcommunity=community_subtype)                               #This generates a data frame by combining all the meatadata I ran the for loop on.


###INITIALIZING TAXONOMY AND PHYLOSEQ OBJECT GENERATION
taxa_table<-read.table("GTDB_RDP.tsv", header=TRUE)                             #This reads the taxonomy table from GTDB by the femmicro pipeline
write.csv(taxa_table,file = "taxa_table.csv")
taxa_table<-as.matrix(taxa_table)                                               #Reading it as a matrix instead of a df
getwd()

rownames(samdf)<-samdf$ID                                                       #Phyloseq won't like if the sample names does not match the seqtab file, so renaming the rownames to the ID fixes the problem of "sample names do not match"
ps<-phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
             sample_data(samdf),tax_table(taxa_table))
dna<-Biostrings::DNAStringSet(taxa_names(ps))
names(dna)<-taxa_names(ps)
ps<-merge_phyloseq(ps,dna)
taxa_names(ps)<-paste0("ASV",seq(ntaxa(ps)))
#taxa_names(ps)

#View(seqtab.nochim)

###GETTING METRICS ON READS
reads_table<-as.data.frame(read.table("Nreads.tsv"))
reads_table$
samdf
total_raw_reads<-sum(as.numeric(reads_table$num_seqs))
total_nochim_reads<-sum(as.numeric(reads_table$nonchim))
total_raw_reads
total_nochim_reads

###COUNTING TAXA AND READS ACROSS PS OBJECTS AND SUBSETS
ps
sum(colSums(ps@otu_table)>1)                                                    #This tells me the # of taxa with more than 1 read per ASV.
length(taxa_names(ps))                                                          #This tells me the number of reference sequences originally there.
#ps_pruned<-prune_taxa(taxa_sums(ps)>0, ps)                                     #I commented this out as all of the ASVs here were present at least once in one sample. 
#ps_pruned

#Subsetting by sex 
ps_female<-subset_samples(ps, Sex=="Female")                                    #This generates a ps object only for the ones that match the female sex. 
ps_male<-subset_samples(ps, Sex=="Male")                                        #This generates a ps object only for the ones that match the male sex.
#Subsetting by community
ps_seroconverted<-subset_samples(ps, Community=="IAB+")                         #I repeat this in every step to stratify it accordingly
ps_nonseroconverted<-subset_samples(ps, Community=="IAB-")
#Subsetting by community AND sex
ps_seroconverted_male<-subset_samples(ps, Sex=="Male" & Community=="IAB+")
ps_seroconverted_female<-subset_samples(ps, Sex=="Female" & Community=="IAB+")
ps_nonseroconverted_male<-subset_samples(ps, Sex=="Male" & Community=="IAB-")
ps_nonseroconverted_female<-subset_samples(ps, Sex=="Female" & Community=="IAB-")
#Subsetting by subcommunity
ps_S2<-subset_samples(ps, Subcommunity=="S2")
ps_S5<-subset_samples(ps, Subcommunity=="S5")
ps_NS1<-subset_samples(ps, Subcommunity=="NS1")
ps_NS6<-subset_samples(ps, Subcommunity=="NS6")
#Subsetting by subcommunity AND sex
ps_S2_male<-subset_samples(ps, Sex=="Male" & Subcommunity=="S2")
ps_S2_female<-subset_samples(ps, Sex=="Female" & Subcommunity=="S2")
ps_S5_male<-subset_samples(ps, Sex=="Male" & Subcommunity=="S5")
ps_S5_female<-subset_samples(ps, Sex=="Female" & Subcommunity=="S5")
ps_NS1_male<-subset_samples(ps, Sex=="Male" & Subcommunity=="NS1")
ps_NS1_female<-subset_samples(ps, Sex=="Female" & Subcommunity=="NS1")
ps_NS6_male<-subset_samples(ps, Sex=="Male" & Subcommunity=="NS6")
ps_NS6_female<-subset_samples(ps, Sex=="Female" & Subcommunity=="NS6")

ps_seroconverted_pruned<-prune_taxa(taxa_sums(ps_seroconverted)>0, ps_seroconverted)                #This creates a new object with the ASVs unique to seroconverted samples, with taxa greater ASV>0 
ps_nonseroconverted_pruned<-prune_taxa(taxa_sums(ps_nonseroconverted)>0, ps_nonseroconverted)       #I repeat this across each sub-stratified ps object.
ps_male_pruned<-prune_taxa(taxa_sums(ps_male)>0, ps_male)
ps_female_pruned<-prune_taxa(taxa_sums(ps_female)>0, ps_female)
ps_S2_pruned<-prune_taxa(taxa_sums(ps_S2)>0, ps_S2)
ps_S5_pruned<-prune_taxa(taxa_sums(ps_S5)>0, ps_S5)
ps_NS1_pruned<-prune_taxa(taxa_sums(ps_NS1)>0, ps_NS1)
ps_NS6_pruned<-prune_taxa(taxa_sums(ps_NS6)>0, ps_NS6)

counting_asvs_per_sample<-function(ps_object){                                                      #This function allows me to get the number of ASVs per ps object, I use it for each object.
  asv_counts <- apply(t(otu_table(ps_object)), 2, function(x) sum(x > 0))
  asv_counts <- data.frame(ASV_counts = asv_counts)
  return(asv_counts)
}

asvs_seroconverted_pruned<-counting_asvs_per_sample(ps_seroconverted_pruned)                        #Running this function across each ps object.
asvs_nonseroconverted_pruned<-counting_asvs_per_sample(ps_nonseroconverted_pruned)
asvs_male_pruned<-counting_asvs_per_sample(ps_male_pruned)
asvs_female_pruned<-counting_asvs_per_sample(ps_female_pruned)
asvs_S2_pruned<-counting_asvs_per_sample(ps_S2_pruned)
asvs_S5_pruned<-counting_asvs_per_sample(ps_S5_pruned)
asvs_NS1_pruned<-counting_asvs_per_sample(ps_NS1_pruned)
asvs_NS6_pruned<-counting_asvs_per_sample(ps_NS6_pruned)

easy_printing<-function(asvs_dataframe){                                        #This prints the output from the function in a nice way so I can copy and paste it in Prism.
  for(item in asvs_dataframe){
    cat(item,"\n")
  }
}

easy_printing(asvs_seroconverted_pruned)                                        #Running the easy_printing function across each object.
easy_printing(asvs_nonseroconverted_pruned)
easy_printing(asvs_male_pruned)
easy_printing(asvs_female_pruned)
easy_printing(asvs_S2_pruned)
easy_printing(asvs_S5_pruned)
easy_printing(asvs_NS1_pruned)
easy_printing(asvs_NS6_pruned)

###ESTIMATING THE NUMBER OF PHYLA, CLASS, ORDER, FAMILY, GENUS, AND SPECIES ACROSS EACH PS OBJECT
taxa_male<-as.data.frame(tax_table(ps_male))
taxa_female<-as.data.frame(tax_table(ps_female))
taxa_male


###ALPHA DIVERSITY ACROSS NON-NORMALIZED DATA
alpha_diversity_male<-estimate_richness(ps_male, measures=c("Shannon","Simpson","Observed","Chao1"))
alpha_diversity_female<-estimate_richness(ps_female, measures=c("Shannon","Simpson","Observed","Chao1"))
alpha_diversity_seroconverted<-estimate_richness(ps_seroconverted, measures=c("Shannon","Simpson","Observed","Chao1"))
alpha_diversity_nonseroconverted<-estimate_richness(ps_nonseroconverted, measures=c("Shannon","Simpson","Observed","Chao1"))
alpha_diversity_S2<-estimate_richness(ps_S2, measures=c("Shannon","Simpson","Observed","Chao1"))
alpha_diversity_S5<-estimate_richness(ps_S5, measures=c("Shannon","Simpson","Observed","Chao1"))
alpha_diversity_NS1<-estimate_richness(ps_NS1, measures=c("Shannon","Simpson","Observed","Chao1"))
alpha_diversity_NS6<-estimate_richness(ps_NS6, measures=c("Shannon","Simpson","Observed","Chao1"))

cat(alpha_diversity_male$Observed,"\n")
cat(alpha_diversity_female$Observed,"\n")
cat(alpha_diversity_seroconverted$Observed,"\n")
cat(alpha_diversity_nonseroconverted$Observed,"\n")
cat(alpha_diversity_S2$Observed,"\n")
cat(alpha_diversity_S5$Observed,"\n")
cat(alpha_diversity_NS1$Observed,"\n")
cat(alpha_diversity_NS6$Observed,"\n")

cat(alpha_diversity_male$Shannon,"\n")
cat(alpha_diversity_female$Shannon,"\n")
cat(alpha_diversity_seroconverted$Shannon,"\n")
cat(alpha_diversity_nonseroconverted$Shannon,"\n")
cat(alpha_diversity_S2$Shannon,"\n")
cat(alpha_diversity_S5$Shannon,"\n")
cat(alpha_diversity_NS1$Shannon,"\n")
cat(alpha_diversity_NS6$Shannon,"\n")

cat(alpha_diversity_male$Simpson,"\n")
cat(alpha_diversity_female$Simpson,"\n")
cat(alpha_diversity_seroconverted$Simpson,"\n")
cat(alpha_diversity_nonseroconverted$Simpson,"\n")
cat(alpha_diversity_S2$Simpson,"\n")
cat(alpha_diversity_S5$Simpson,"\n")
cat(alpha_diversity_NS1$Simpson,"\n")
cat(alpha_diversity_NS6$Simpson,"\n")


#closeAllConnections()
###NORMALIZING THE PS OBJECTS AND GETTING THE BETA-DIVERSITY PLOT
ps_clr<-norm_clr(ps)

graphing_beta_diversity<-function(ps_all_samples_minus_empty,method_plot,dist,string,category){
  getwd()
  sink(paste0(method_plot,"_",string,"_",dist,"_",category,".txt"))
  group_colours<-c("skyblue","mediumpurple")
  beta_dist_ord <- ordinate(ps_all_samples_minus_empty, method = method_plot, distance = dist)
  beta_dist_plot <- plot_ordination(ps_all_samples_minus_empty, ordination = beta_dist_ord, type = "ID", color = "Sex")
  print(beta_dist_ord)
  plot <- beta_dist_plot + scale_color_manual(values=group_colours) + 
    stat_ellipse(alpha = 0.07, geom = "polygon", aes(fill = Sex),show.legend=FALSE) +
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
  print(adonis2(distance_measure ~ Sex, data=sampledf, method=dist))
  print(plot)
  sink()
  ggsave(paste0(method_plot,"_",dist,"_",string,"_",category,".pdf"),width=5.5,height=3.5)

}

sample_data(ps)
graphing_beta_diversity(ps_clr,"NMDS","bray","Week10","sex")
graphing_beta_diversity(ps_clr,"NMDS","jaccard","Week10","sex")
graphing_beta_diversity(ps_clr,"PCoA","bray","Week10","sex")
graphing_beta_diversity(ps_clr,"PCoA","jaccard","Week10","sex")
graphing_beta_diversity(ps_clr,"PCoA","euclidean","Week10","sex")

graphing_beta_diversity2<-function(ps_all_samples_minus_empty,method_plot,dist,string,category){
  getwd()
  sink(paste0(method_plot,"_",string,"_",dist,"_",category,".txt"))
  group_colours<-c("green4","blue")
  beta_dist_ord <- ordinate(ps_all_samples_minus_empty, method = method_plot, distance = dist)
  beta_dist_plot <- plot_ordination(ps_all_samples_minus_empty, ordination = beta_dist_ord, type = "ID", color = "Community")
  print(beta_dist_ord)
  plot <- beta_dist_plot + scale_color_manual(values=group_colours) + 
    stat_ellipse(alpha = 0.07, geom = "polygon", aes(fill = Community),show.legend=FALSE) +
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
  print(adonis2(distance_measure ~ Community, data=sampledf, method=dist))
  print(plot)
  sink()
  ggsave(paste0(method_plot,"_",dist,"_",string,"_",category,".pdf"),width=5.5,height=3.5)
}

graphing_beta_diversity2(ps_clr,"NMDS","bray","Week10","community")
graphing_beta_diversity2(ps_clr,"NMDS","jaccard","Week10","community")
graphing_beta_diversity2(ps_clr,"PCoA","bray","Week10","community")
graphing_beta_diversity2(ps_clr,"PCoA","jaccard","Week10","community")
graphing_beta_diversity2(ps_clr,"PCoA","euclidean","Week10","community")

graphing_beta_diversity3<-function(ps_all_samples_minus_empty,method_plot,dist,string,category){
  getwd()
  sink(paste0(method_plot,"_",string,"_",dist,"_",category,".txt"))
  group_colours<-c("gray60","cadetblue3","dodgerblue3","khaki3")
  beta_dist_ord <- ordinate(ps_all_samples_minus_empty, method = method_plot, distance = dist)
  beta_dist_plot <- plot_ordination(ps_all_samples_minus_empty, ordination = beta_dist_ord, type = "ID", color = "Subcommunity")
  print(beta_dist_ord)
  plot <- beta_dist_plot + scale_color_manual(values=group_colours) + 
    stat_ellipse(alpha = 0.07, geom = "polygon", aes(fill = Subcommunity),show.legend=FALSE) +
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
  print(pairwise.adonis2(distance_measure ~ Subcommunity, data=sampledf, method=dist,nperm = 999))
  print(plot)
  sink()
  ggsave(paste0(method_plot,"_",dist,"_",string,"_",category,".pdf"),width=5.5,height=3.5)
}

graphing_beta_diversity3(ps_clr,"NMDS","bray","Week10","subcommunity")
graphing_beta_diversity3(ps_clr,"NMDS","jaccard","Week10","subcommunity")
graphing_beta_diversity3(ps_clr,"PCoA","bray","Week10","subcommunity")
graphing_beta_diversity3(ps_clr,"PCoA","jaccard","Week10","subcommunity")
graphing_beta_diversity3(ps_clr,"PCoA","euclidean","Week10","subcommunity")

###ENGRAFTMENT ANALYSES, WHICH SPECIES ENGRAFTED BY THE WEEKS THEY ARE ANALYZED?
unique_phyla_seroconverted<-get_taxa_unique(ps_seroconverted_pruned, taxonomic.rank = rank_names(ps)[2],errorIfNULL = TRUE)
unique_phyla_seroconverted
unique_phyla_nonseroconverted<-get_taxa_unique(ps_nonseroconverted_pruned, taxonomic.rank = rank_names(ps)[2],errorIfNULL = TRUE)
unique_phyla_nonseroconverted
unique_phyla_NS1<-get_taxa_unique(ps_NS1_pruned, taxonomic.rank = rank_names(ps)[2],errorIfNULL = TRUE)
unique_phyla_NS1
unique_phyla_NS6<-get_taxa_unique(ps_NS6_pruned, taxonomic.rank = rank_names(ps)[2],errorIfNULL = TRUE)
unique_phyla_NS6
unique_phyla_S2<-get_taxa_unique(ps_S2_pruned, taxonomic.rank = rank_names(ps)[2],errorIfNULL = TRUE)
unique_phyla_S2
unique_phyla_S5<-get_taxa_unique(ps_S5_pruned, taxonomic.rank = rank_names(ps)[2],errorIfNULL = TRUE)
unique_phyla_S5


###PAINTER PLOTS. BY PHYLA AND BY GENUS
abundance_plots<-function(ps_object, top_to_look_for, type_of_plot, norm.method){
  clrs <- c("#e3869c","#5db746","#a25acd","#a3b334","#5d6bc8","#d29f36","#cf4ea6","#5bbe7d","#d4426c",
            "#4bb8aa","#cf4336","#5f96d0","#d96e2d","#bd8bd4","#487c3d","#974b7e","#a9af62","#a85147",
            "#7e6e28","#cc8957")
  file_name<-paste0(type_of_plot,"_",top_to_look_for,"_",norm.method)
  ps_object.genus<-tax_glom(ps_object, taxrank = "Phylum")
  top<-names(sort(taxa_sums(ps_object.genus), decreasing=TRUE))[1:top_to_look_for]
  if(type_of_plot=="absolute"){
    ps_absolute_abundance<-prune_taxa(top, ps_object.genus)
    plot_genus_top<-plot_bar(ps_absolute_abundance, fill="Phylum") + scale_fill_manual(values=clrs)
    plot_genus_top + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity") + scale_color_manual(values=clrs) + theme(
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
    plot_genus_top <- plot_bar(ps_relative_abundance.genus, fill = "Phylum") +
      scale_fill_manual(values = clrs) +
      geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity") +
      scale_color_manual(values = clrs) +
      theme(
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "gray", size = 0.20),
        panel.grid.minor = element_line(colour = "gray", size = 0.05),
      )
    otu_df<-as.data.frame(otu_table(ps_relative_abundance.genus))
    otu_df<-otu_df %>%
      mutate(total_relative_abundance = rowSums(.))
    otu_df<- otu_df %>%
      mutate(new_column = 100 - total_relative_abundance)
    write.csv(x = otu_df, file=paste0(file_name,"_ASVcount.csv"))
    write.csv(x = tax_table(ps_relative_abundance.genus), file=paste0(file_name,"_taxaNames.csv"))
  }
  ggsave(filename = paste0(file_name,".pdf"),plot = plot_genus_top)
}
getwd()
abundance_plots(ps_seroconverted_pruned,7,"absolute","seroconverted_w9w10")
abundance_plots(ps_seroconverted_pruned,7,"relative","seroconverted_w9w10")
abundance_plots(ps_nonseroconverted_pruned,7,"absolute","nonseroconverted_w9w10")
abundance_plots(ps_nonseroconverted_pruned,7,"relative","nonseroconverted_w9w10")
abundance_plots(ps_female_pruned,7,"relative","female_w9w10")
abundance_plots(ps_male_pruned,7,"relative","male_w9w10")
abundance_plots(ps_male_pruned,7,"absolute","male_w9w10")
abundance_plots(ps_female_pruned,7,"absolute","female_w9w10")

View(tax_table(ps_nonseroconverted_pruned))


taxa_names()#Differential abundance of metagenome prediction
maaslin_function<-function(ps_object,type_to_analyze,fixed_value){
  metadata<-as.data.frame(sample_data(ps_object))
  class(metadata)<-"data.frame"
  ps_otu<-as(otu_table(ps_object), "matrix")
  ps_otu_df<-as.data.frame(ps_otu)
  fit_data = Maaslin2(input_data = ps_otu_df, 
                      input_metadata = metadata, 
                      output=type_to_analyze,
                      fixed_effects=fixed_value,
                      normalization = "NONE")                                   #Change this to the variable you want to make a comparison with
  print(fit_data$results)
}
sample_data(ps_clr)
maaslin_function(ps_clr,"maaslin_da_sex","Sex")
maaslin_function(ps_clr,"maaslin_da_community","Community")



###WORKING ON PICRUST2 DATA AND MAASLIN2
#This one is cross-referenced to the default run from Picrust2 (not local mode)
picrust_to_maaslin2<-function(path_picrust2, file_to_use, variable, normalization_method, transformation){
  setwd(path_picrust2)
  #maaslin_pathways_df<-(read.delim("pathways_out/path_abun_unstrat_renamed.txt",header = TRUE))                                  #I had to rename the file, as the pathways names were generating an error in Maaslin2
  maaslin_pathways_df<-read.delim(file_to_use)
  colnames(maaslin_pathways_df) <- sub("^.", "", colnames(maaslin_pathways_df), perl = TRUE)                                 #R was adding an 'X' at the beginning of my sample names, so I removed that 'X'.
  maaslin_pathways_df<-t(maaslin_pathways_df)
  colnames(maaslin_pathways_df)<-maaslin_pathways_df[1,]
  maaslin_pathways_df<-maaslin_pathways_df[-1,]
  maaslin_pathways_df<-as.data.frame(maaslin_pathways_df)
  maaslin_pathways_df[] <- lapply(maaslin_pathways_df, as.numeric)
  metadata<-as.data.frame(sample_data(ps))
  class(metadata)<-"data.frame"
  write.matrix(x=metadata,file = "metadata_file.txt")
  fit_data = Maaslin2(input_data = maaslin_pathways_df, input_metadata = metadata, output="maaslin2_results_default_tss_log2", 
                    fixed_effects=variable, normalization=normalization_method,transform = transformation)
}


picrust_to_maaslin2(path_picrust2 = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2/picrust_july2024/default_run/default_output/pathways_out/",
                    variable = "Community", normalization_method = "TSS", transformation = "LOG", file_to_use ="path_abun_unstrat_renamed.txt" )               #Make sure to rename the files if needed
picrust_to_maaslin2(path_picrust2 = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2/picrust_july2024/local_output_run/EC_metagenome_out/",
                    variable = "Community", normalization_method = "TSS", transformation = "LOG")
picrust_to_maaslin2(path_picrust2 = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2/picrust_july2024/default_run/default_output/KO_metagenome_out/",
                    variable = "Community", normalization_method = "TSS", transformation = "LOG", file_to_use = "pred_metagenome_unstrat.tsv")
picrust_to_maaslin2(path_picrust2 = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2/picrust_july2024/local_output_run/KO_metagenome_out/",
                    variable = "Community", normalization_method = "TSS", transformation = "LOG", file_to_use = "pred_metagenome_unstrat.tsv")
picrust_to_maaslin2(path_picrust2 = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2/picrust_july2024/local_output_run/KO_metagenome_out/",
                    variable = "Community", normalization_method = "TSS", transformation = "LOG", file_to_use = "pred_metagenome_unstrat.tsv")

###HEATMAP VISUALIZATION OF THE COVERAGE OF EACH PATHWAY FROM METACYC (LINKED WITH PICRUST2)
heatmap_function<-function(file_to_use){
  coverage_pathway_df<-read.table(file_to_use,header = TRUE)
  rownames(coverage_pathway_df)<-coverage_pathway_df[,1]
  coverage_pathway_df<-coverage_pathway_df[,-1]
  coverage_pathway_df<-as.matrix(coverage_pathway_df)
  coverage_pathway_df
  ncol(coverage_pathway_df)
  my_palette <- colorRampPalette(c("white", "blue"))(n = 300)
  heatmap.2(coverage_pathway_df,
        col = my_palette,
        Colv = "none",
        trace = "none",
        density.info="none",
        ColSideColors = c(rep("#0077FF",68), rep("green",101)))
}
heatmap_function("../../default_run/default_output/pathway_coverage_unstrat_sorted.txt")
heatmap_function("../../local_output_run/EC_metagenome_out/maaslin2_results_default_tss_log2/path_cov_unstrat_sorted.txt")

###WORKING ON MICROBIOME DATA, PICRUST2 DATA AND DESEQ2 DIFFERENTIAL ABUNDANCE
sample_data(ps)
head(sample_data(ps)$Community)                                                 #Sanity check, all are IAB+ or IAB-
deseq2_object = phyloseq_to_deseq2(ps, ~ Sex)
deseq2_object = DESeq2::DESeq(deseq2_object, test="Wald",fitType = "mean")
