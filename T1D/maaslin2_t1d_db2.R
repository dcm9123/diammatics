#Maaslin2
#July 21st, 2025

library(Maaslin2)

picrust_to_maaslin2<-function(path_picrust2, file_to_use, variable, normalization_method,maaslin_path, 
                              transformation,metadata_file,feature_to_use,random_input,reference_val){
  setwd(path_picrust2)
  maaslin_input = read.table(file_to_use, header = TRUE, row.names = 1, sep='\t')
  metadata_input = read.csv(metadata_file, header = TRUE, row.names = 1)
  output_file = paste0(maaslin_path,"maaslin2_",normalization_method,"_",transformation,"_",variable,"_",feature_to_use)
  fixed_effects_to_use = variable
  transformation_method = transformation
  normalization_to_use = normalization_method
  random_to_use = random_input
  reference_value = reference_val
  fit_data = Maaslin2(input_data = maaslin_input, input_metadata = metadata_input, output=output_file,
                      fixed_effects = fixed_effects_to_use, transform = transformation_method, 
                      normalization = normalization_to_use, random_effects = random_to_use, reference = reference_value)
  
  cat("The results were generated based on the following information: \n")
  cat(paste0("Input file: ",path_picrust2,file_to_use,"\n"))
  cat(paste0("Metadata file: ",metadata_file,"\n"))
  cat(paste0("Fixed effect variable: ",variable,"\n"))
  cat(paste0("Random effect variable: ",random_input,"\n"))
  cat(paste0("Reference value: ",reference_val,"\n"))
  cat(paste0("Feature: ",feature_to_use,"\n"))
  cat(paste0("Normalization method: ",normalization_method,"\n"))
  cat(paste0("Transformation method: ",transformation,"\n"))
  cat(paste0("File written to: ",maaslin_path,output_file,"\n"))
  cat("\n")
  cat("Finished!")
}


#Comparing week and consortia (NS1, NS6, S2, S5) and w5w6, w9w10. Reference NS1_w5w6

#Metadata variables: ID,Sex,Community,Subcommunity,Timepoint,Merged_weeks,Week_and_consortia
picrust_to_maaslin2(path_picrust2 = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2_june232025/",
                    variable = "Sex", 
                    normalization_method = "TSS", 
                    transformation = "LOG", 
                    file_to_use ="KO_merged_metagenome.tsv",
                    metadata_file = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv",
                    feature_to_use = "KO",
                    random_input = "",
                    reference_val = "Sex,Male", 
                    maaslin_path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/") 

picrust_to_maaslin2(path_picrust2 = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2_june232025/",
                    variable = "Merged_weeks", 
                    normalization_method = "TSS", 
                    transformation = "LOG", 
                    file_to_use ="KO_merged_metagenome.tsv",
                    metadata_file = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv",
                    feature_to_use = "KO",
                    random_input = "Sex",
                    reference_val = "",
                    maaslin_path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/")

picrust_to_maaslin2(path_picrust2 = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2_june232025/",
                    variable = "Subcommunity", 
                    normalization_method = "TSS", 
                    transformation = "LOG", 
                    file_to_use ="KO_merged_metagenome.tsv",
                    metadata_file = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv",
                    feature_to_use = "KO",
                    random_input = "Sex",
                    reference_val = "Subcommunity,NS1",
                    maaslin_path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/")

picrust_to_maaslin2(path_picrust2 = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2_june232025/",
                    variable = "Community", 
                    normalization_method = "TSS", 
                    transformation = "LOG", 
                    file_to_use ="KO_merged_metagenome.tsv",
                    metadata_file = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv",
                    feature_to_use = "KO",
                    random_input = "Sex",
                    reference_val = "",
                    maaslin_path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/")

picrust_to_maaslin2(path_picrust2 = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2_june232025/",
                    variable = "Week_and_consortia", 
                    normalization_method = "TSS", 
                    transformation = "LOG", 
                    file_to_use ="EC_merged_metagenome.tsv",
                    metadata_file = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv",
                    feature_to_use = "EC",
                    random_input = "Sex",
                    reference_val = "Week_and_consortia,NS1_w5w6", 
                    maaslin_path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/")  

picrust_to_maaslin2(path_picrust2 = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2_june232025/",
                    variable = "Week_and_consortia", 
                    normalization_method = "TSS", 
                    transformation = "LOG", 
                    file_to_use ="KO_merged_metagenome.tsv",
                    metadata_file = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv",
                    feature_to_use = "KO",
                    random_input = "Sex",
                    reference_val = "Week_and_consortia,NS1_w5w6", 
                    maaslin_path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/") 

#Comparing week and consortia using S2_w5w6 as reference
picrust_to_maaslin2(path_picrust2 = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2_june232025/",
                    variable = "Week_and_consortia", 
                    normalization_method = "TSS", 
                    transformation = "LOG", 
                    file_to_use ="KO_merged_metagenome.tsv",
                    metadata_file = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv",
                    feature_to_use = "KO",
                    random_input = "Sex",
                    reference_val = "Week_and_consortia,S2_w5w6",
                    maaslin_path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/")

picrust_to_maaslin2(path_picrust2 = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2_june232025/",
                    variable = "Week_and_consortia", 
                    normalization_method = "TSS", 
                    transformation = "LOG", 
                    file_to_use ="EC_merged_metagenome.tsv",
                    metadata_file = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv",
                    feature_to_use = "EC",
                    random_input = "Sex",
                    reference_val = "Week_and_consortia,S2_w5w6",
                    maaslin_path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/")
