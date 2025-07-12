# PICRUSt2
For PICRUST2 I have to create 4 individual phylogenetic trees from the consortia. For this, I will use the files located in `/bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/16S_data/16S_after_fixing_and_discussion/good_files_to_use/t1d_consortia_dbs` , specifically the ones called:

`rm_dupt1d_db2_withDups_S_NS1_species.fasta`
`rm_dupt1d_db2_withDups_S_NS6_species.fasta`
`rm_dupt1d_db2_withDups_S_S2_species.fasta`
`rm_dupt1d_db2_withDups_S_S5_species.fasta`

## MAFFT
These files have the duplicate sequences removed of EACH consortia, not from the global t1d updated db 2. So, the only duplicates that were removed where located in each consortia. 

The next step is to run MAFFT in each file:
`mafft --thread 28 --auto --reorder "t1d_consortia_dbs/rm_dupt1d_db2_withDups_S_S2_species.fasta" > "MAFFT_consortia/rm_dupt1d_db2_withDups_S_S2_species_aligned.fasta"`

`mafft --thread 28 --auto --reorder "t1d_consortia_dbs/rm_dupt1d_db2_withDups_S_S5_species.fasta" > "MAFFT_consortia/rm_dupt1d_db2_withDups_S_S5_species_aligned.fasta"`

`mafft --thread 28 --auto --reorder "t1d_consortia_dbs/rm_dupt1d_db2_withDups_S_NS6_species.fasta" > "MAFFT_consortia/rm_dupt1d_db2_withDups_S_NS6_species_aligned.fasta"`

`mafft --thread 28 --auto --reorder "t1d_consortia_dbs/rm_dupt1d_db2_withDups_S_NS1_species.fasta" > "MAFFT_consortia/rm_dupt1d_db2_withDups_S_NS1_species_aligned.fasta"`

Make sure to remove identical header names and identical sequences within each fasta file! RAxML does not process the tree if duplicates are found!

Each one of the aligned files can be found in `/bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/16S_data/16S_after_fixing_and_discussion/good_files_to_use/MAFFT_consortia`

### June 26th, 2025
Yesterday I was told about a bug in the latest FemMicro16S pipeline (the ones I have in the blackbox folder in ARC), where one of the taxa had an old Phyla nomenclature - Firmicutes. I fixed this problem by assigning the right phylum name, and I also made sure the 082 and 086 sequences were properly mapping to the taxa pointed by FemMicro. I still have to review the files where the duplicates have been removed, as we are currently running FemMicro with the duplicates included (duplicated sequences). Then I ran the seqtabnochim objects into my R code (which I just updated) so I could extract and merge the 5 ps objects. We obtained 535 taxa from this. I need to check, however, the sequence that I merge manually in R into the taxa file, which does not match the ASV_ID_Seq.

In addition to that, I extracted the latest metrics from the pipeline and shared them with Laura. Some of these metrics are:

`/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/t1d_db_fixed_discussed/FemMicro_Daniel/Nreads_5plates.tsv`
`/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/t1d_db_fixed_discussed/FemMicro_Daniel/cutAdapt_metrics_5plates.xlsx`
`/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/t1d_db_fixed_discussed/FemMicro_Daniel/plate1.1/multiqc_report_filtered_p1.html` (Repeated for each plate, giving me a total of 5).

## RAxML

For RAxML, I had this code that gets bootstrap trees, followed by the selection of the best tree. The code is called `raxml_slurm.sh`, and it is found in this repository and in ARC at `/bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/RAxML
` 
In this part, I had to make sure that the headers of the sequences didn't have any spaces, and to assign just 1 or 2 threads to the job, as more than that will crash the code. Specifically in the `raxml-ng --search` part, I had to include `--force` to force it to run with 1 thread.

### July 1st, 2025
Ignore what I previously wrote. Turns out I was not specifying the number of threads properly for the `raxml-ng --search` part, and when that happens, the program will try to use whatever number of threads the node has available. Look at the updated raxml_slurm.sh file I have added to this repo in this folder in case you need it again.

### July 3rd, 2025
A total of 17 genomes were flagged by NCBI. 14 of those had illumina adapters in one to three contigs in each file. Each sequence that needed to be removed was shorter than 60 bp. Laura and I discussed and we decided it's best to manually remove these short sequences, reannotate these genomes, and re-run the quality metrics (checkM, quast, etc). For the remaining 4 genomes, these 4 have contaminants from other species. We removed the contigs associated with these others and I'll be running the same metrics as the other genomes. The command used for Prokka is found in this folder.

### July 8th, 2025

After a few days of coding I was able to come up with a finalized script that takes the raw eggnog_annotations file after running eggnog on various samples, and then it produces two final files containing the KOs and ECs formatted for PICRUSt2 to read properly. For this to happen, the user must have and provide:

1) The number of path(s) where their samples are located, as well as the path(s) itself.
2) The formatting of the file name must follow this nomenclature: `S_NS1_Am_001.emapper.annotations`. The first letter indicating it's a sample ("S"), the second instance indicating the consortia name/identifier (`NS1 = Nonseroconverted 1`), the third instance a microbe identifier (in this case, `Am = Akkermansia muciniphila`), and fourth, a unique number identifier (`001 = sample # 001)`. It must also end in `emapper.annotations` for it to be read.
3) Provide an output path to store the output files.
4) Have the following packages installed:
  `pandas (v1.4.4)`,
  `tqdm (v4.65.0)`,
  `python (v3.11.5)`

The name of the script is: `eggnog_to_picrust2.py`, and the two examples final files it generates are `final_EC_file.tsv` and `Final_KO_file.tsv`

###July 9th, 2025

I have downloaded the PICRUSt2 original files so I can used them as a reference format for mine. The files are called: `16S_reference.txt`, `ec_reference.txt.gz`, and `ko_reference_cut.txt` . Some of these files were too large to be uploaded, so I cut the ko_reference to the top 5 lines, and I left the EC file gzipped.
I need to add the word 'assembly' instead of Sample ID to the final formatted EC and KO files from my `eggnog_to_picrust2.py` script. I also need to leave the 'ko:' in the KO identifiers so it looks like `'ko:K02237'` .On the other hand, the 'EC' word should be removed, in a way that it only displays the number (i.e. `2.3.5.16`).
The 16S.txt file must have the word `assembly` for the genome identifiers, a tab, and the column name 16S_rRNA_count as well.

###July 11th, 2025
I have restarted the Anvi'o analysis for the new genomes. For this, I had to start over using Anvi'o again. Some of the things that I've learned whilst installing the anvio-dev version, is that the databases are not saved in a specific directory of your choice (or if there is, it's more complicated than I thought...). Instead, it is saved in the same directory where Anaconda is found (in my case it was `/home/daniel.castanedamogo/github`). Then, I had to do the following:

1. Identify the paths where my genomes are located, and make a new file with the sample id and the path where it is located: `/bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/anvio/all_genomes_NCBI/name_and_path_all.txt`.
2. Download the default config file by running `anvi-run-workflow -w contigs --get-default-config config-contigs-default.json`
3. Download the databases I am planning to use (in my case KEGG, COG, and Taxonomy). For this, I ran: `anvi-setup-ncbi-cogs --threads 16`, followed by `anvi-setup-scg-taxonomy` (this one downloads GTDB genome information), and finally `anvi-setup-kegg-data`. The databases will be downloaded in the home directory, so make sure you have at least ~200-300 GB free.
4. Modify the config file, in a way that it sets the number of threads, and the file containing the path and name of the samples is named properly.
5. For some reason, the contigs worklow kept crashing when submitting it through slurm, so I decided to use an interactive node instead.
6. Run the contigs workflow like this: `anvi-run-workflow -w contigs -c config_default.json`. In my case, with 116 samples, I had to run the interactive node 3 times for 5 hours each. However, when rerunning it, anvi'o is a bit picky on rewriting some of the files and pick up where it left off. So, in that case I had to run this command `anvi-run-workflow - contigs -c config_default.json --additional-params --rerun-incomplete --keep-going`. This gave me the contigs database I need for downstream steps.



