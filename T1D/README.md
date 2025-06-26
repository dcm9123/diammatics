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

