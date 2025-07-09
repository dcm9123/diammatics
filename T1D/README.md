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
  `pandas (v1.4.4)
  tqdm (v4.65.0)
  python (v3.11.5)`

