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

### July 9th, 2025

### **<ins> Anvi'o </ins>**

I have downloaded the PICRUSt2 original files so I can used them as a reference format for mine. The files are called: `16S_reference.txt`, `ec_reference.txt.gz`, and `ko_reference_cut.txt` . Some of these files were too large to be uploaded, so I cut the ko_reference to the top 5 lines, and I left the EC file gzipped.
I need to add the word 'assembly' instead of Sample ID to the final formatted EC and KO files from my `eggnog_to_picrust2.py` script. I also need to leave the 'ko:' in the KO identifiers so it looks like `'ko:K02237'` .On the other hand, the 'EC' word should be removed, in a way that it only displays the number (i.e. `2.3.5.16`).
The 16S.txt file must have the word `assembly` for the genome identifiers, a tab, and the column name 16S_rRNA_count as well.

### July 11th, 2025

### **<ins> Anvi'o </ins>**

I have restarted the Anvi'o analysis for the new genomes. For this, I had to start over using Anvi'o again. Some of the things that I've learned whilst installing the anvio-dev version, is that the databases are not saved in a specific directory of your choice (or if there is, it's more complicated than I thought...). Instead, it is saved in the same directory where Anaconda is found (in my case it was `/home/daniel.castanedamogo/github`). Then, I had to do the following:

1. Identify the paths where my genomes are located, and make a new file with the sample id and the path where it is located: `/bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/anvio/all_genomes_NCBI/name_and_path_all.txt`.
2. Download the default config file by running `anvi-run-workflow -w contigs --get-default-config config-contigs-default.json`
3. Download the databases I am planning to use (in my case KEGG, COG, and Taxonomy). For this, I ran: `anvi-setup-ncbi-cogs --threads 16`, followed by `anvi-setup-scg-taxonomy` (this one downloads GTDB genome information), and finally `anvi-setup-kegg-data`. The databases will be downloaded in the home directory, so make sure you have at least ~200-300 GB free.
4. Modify the config file, in a way that it sets the number of threads, and the file containing the path and name of the samples is named properly.
5. For some reason, the contigs worklow kept crashing when submitting it through slurm, so I decided to use an interactive node instead.
6. Run the contigs workflow like this: `anvi-run-workflow -w contigs -c config_default.json`. In my case, with 116 samples, I had to run the interactive node 3 times for 5 hours each. However, when rerunning it, anvi'o is a bit picky on rewriting some of the files and pick up where it left off. So, in that case I had to run this command `anvi-run-workflow -w contigs -c config_default.json --additional-params --rerun-incomplete --keep-going`. This gave me the contigs database I need for downstream steps.

The next part consisted in getting the sequences for hmm hits. For this to happen, I needed to run this short bash script across the new folder generated by the contigs workflow (`02_CONTIGS/`):

`for file in *-contigs.db; do anvi-get-sequences-for-hmm-hits -c ${file} --hmm-source Bacteria_71 -o ${file%%-contigs.db}_genes.fasta; done;`

The next step is to generate a .txt file that describes the hmm sources for a collection of genomes, instead of having individual files for each one of them. For this, I ran:

`anvi-script-gen-hmm-hits-matrix-across-genomes -e name_and_path_to_contigs.txt --hmm-source Bacteria_71 -o hmm_hits_matrix_across_116_genomes.txt`  The output file from this command shows me which genomes have which hmm source.

The next part was to create a new artifact containing the genome storage information that I have so far. For this to be done, I ran: `anvi-gen-genomes-storage -e name_and_path_to_contigs.txt --output-file storage_116-GENOMES.db`. Whilst this is typically ran in a pangenomics workflow, I will be using these genes to generate a gene concatenated file for a phylogenomics tree.

Turns out I did not need the genome storage for the concatenated genes! Instead I had to play around with the parameters until it worked out. However, when getting an error regarding 'file does not exist at `/home/daniel.castanedamogo/tmp/98298asz` make sure to run this command to remove any temporary files related to previous runs:

`find /home/daniel.castanedamogo/tmp -type d -name "tmp*" -exec rm -rf {} +`

Then run this:
```
for file in *-contigs.db; do anvi-get-sequences-for-hmm-hits -c ${file} --hmm-source Bacteria_71 --gene-names ADK,AICARFT_IMPCHas,ATP-synt,ATP-synt_A,Chorismate_synt,EF_TS,Exonuc_VII_L,GrpE,Ham1p_like,IPPT,OSCP,PGK,Pept_tRNA_hydro,RBFA,RNA_pol_L,RNA_pol_Rpb6,RRF,RecO_C,Ribonuclease_P,Ribosom_S12_S23,Ribosomal_L1,Ribosomal_L13,Ribosomal_L14,Ribosomal_L16,Ribosomal_L17,Ribosomal_L18p,Ribosomal_L19,Ribosomal_L2,Ribosomal_L20,Ribosomal_L21p,Ribosomal_L22,Ribosomal_L23,Ribosomal_L27,Ribosomal_L27A,Ribosomal_L28,Ribosomal_L29,Ribosomal_L3,Ribosomal_L32p,Ribosomal_L35p,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6,Ribosomal_L9_C,Ribosomal_S10,Ribosomal_S11,Ribosomal_S13,Ribosomal_S15,Ribosomal_S16,Ribosomal_S17,Ribosomal_S19,Ribosomal_S2,Ribosomal_S20p,Ribosomal_S3_C,Ribosomal_S6,Ribosomal_S7,Ribosomal_S8,Ribosomal_S9,RsfS,RuvX,SecE,SecG,SecY,SmpB,TsaE,UPF0054,YajC,eIF-1a,ribosomal_L24,tRNA-synt_1d,tRNA_m1G_MT,Adenylsucc_synt --get-aa-sequences --concatenate-genes --return-best-hit -o ${file%%-contigs.db}_concatenated_aa.fasta --just-do-it; done;
```

This approach generated a single concatenated fasta file for each sample using the 71 single copy core bacterial genes. However, I do need to get one single concatenated aligned file from all samples in order to create the phylogenomics tree. Everything I've run from getting the sequences from hmm hits until this was not necessary to generate a phylogenomics tree. However, these files could be used for alternative analyses, so it's not so bad to have them!

To gather one single final merged concatenated amino acid file, I ran the following command:
```
anvi-get-sequences-for-hmm-hits --external-genomes ../name_and_path_to_contigs.txt -o final_concatenated_proteins.fa --hmm-source Bacteria_71 --return-best-hit --get-aa-sequences --concatenate --gene-names ADK,AICARFT_IMPCHas,ATP-synt,ATP-synt_A,Chorismate_synt,EF_TS,Exonuc_VII_L,GrpE,Ham1p_like,IPPT,OSCP,PGK,Pept_tRNA_hydro,RBFA,RNA_pol_L,RNA_pol_Rpb6,RRF,RecO_C,Ribonuclease_P,Ribosom_S12_S23,Ribosomal_L1,Ribosomal_L13,Ribosomal_L14,Ribosomal_L16,Ribosomal_L17,Ribosomal_L18p,Ribosomal_L19,Ribosomal_L2,Ribosomal_L20,Ribosomal_L21p,Ribosomal_L22,Ribosomal_L23,Ribosomal_L27,Ribosomal_L27A,Ribosomal_L28,Ribosomal_L29,Ribosomal_L3,Ribosomal_L32p,Ribosomal_L35p,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6,Ribosomal_L9_C,Ribosomal_S10,Ribosomal_S11,Ribosomal_S13,Ribosomal_S15,Ribosomal_S16,Ribosomal_S17,Ribosomal_S19,Ribosomal_S2,Ribosomal_S20p,Ribosomal_S3_C,Ribosomal_S6,Ribosomal_S7,Ribosomal_S8,Ribosomal_S9,RsfS,RuvX,SecE,SecG,SecY,SmpB,TsaE,UPF0054,YajC,eIF-1a,ribosomal_L24,tRNA-synt_1d,tRNA_m1G_MT,Adenylsucc_synt
```

That generated my file, and it showed the parameters and information that was obtaied from the tree:

**Input alignment file path** .....................: `/bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/anvio/all_genomes_NCBI/final_concatenated_proteins.fa` **Output file path** .............................: `/bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/anvio/all_genomes_NCBI/NCBI_phylogenomics_tree.txt` <details> <summary><strong>Alignment names</strong> (click to expand)</summary> S_S5_Rf_112, S_S5_Cc_098, S_S5_Sc_120, S_NS6_Hh_054, S_NS1_Fp_020_S27, S_S5_Bv_093, S_S5_Cb_096, S_S5_Pd_110, S_S5_Cc_097, S_S2_Blss_067_S30, S_S2_Cc_069, S_S5_Ef_104, S_NS6_Fp_052, S_NS6_El_044, S_NS6_Ef_050, S_S2_Fp_079, S_NS1_Bs_004, S_NS1_Bv_007, S_NS1_Ss_026, S_NS1_Fp_019, S_S5_Rg_113, S_NS1_Bu_006, S_NS1_Ao_003, S_S2_Bs_063, S_S2_Bv_065, S_NS1_Ca_028, S_NS1_Shsn_025, S_NS1_Bp_010, S_NS1_Em_018, S_NS6_Ef_049, S_S5_Em_101, S_S2_Et_074_S32, S_NS1_Pd_022, S_S5_Dg_099, S_NS6_Em_046, S_S5_Et_102, S_NS1_El_014, S_S2_Ef_078_S34, S_NS6_Bs_033, S_NS1_Sw_027, S_NS6_Bs_032, S_NS6_Blss_040, S_NS1_Ig_115, S_NS1_Am_001, S_NS6_El_045, S_NS6_Ao_031, S_S2_Ea_077, S_S5_Am_089, S_S5_Pd_111, S_NS1_Cb_011, S_NS1_Et_015, S_NS1_He_021, S_S5_El_100, S_NS6_Bu_117, S_S5_Ea_103, S_S2_He_082_S35, S_S2_Bd_066, S_S2_Ca_088_S36, S_S5_Af_090, S_NS1_Eh_016, S_S2_Rf_084, S_S2_Rg_085, S_NS1_Af_002, S_NS6_Pp_057_S29, S_NS1_Pm_023, S_NS1_Bl_116, S_NS6_Bu_036, S_S2_Am_061, S_S2_Et_075_S33, S_NS6_Bx_038_S28, S_S5_Af_091, S_NS6_He_053, S_NS6_Sw_060, S_NS6_Et_048, S_NS6_Ss_059, S_S2_Bd_064, S_S2_Cs_071, S_S2_Rf_118, S_S2_Fp_080, S_NS6_Cb_041, S_NS1_Ba_008, S_NS6_Bu_034, S_NS1_Cc_012, S_S2_Pd_083, S_S5_Blss_094, S_NS1_Cp_013, S_S5_He_108, S_NS6_Bv_037, S_S5_Lp_109_S38, S_S5_Bp_095, S_NS1_Em_017, S_S2_Ea_076, S_S5_Fp_105, S_NS6_Ba_039, S_NS6_Af_030, S_NS6_Sh_058, S_NS6_Ef_051, S_S5_Fp_106, S_S2_Hf_081, S_NS6_Em_047, S_NS6_Pd_055, S_NS6_Cp_043, S_NS1_Bt_005, S_S2_Sc_086, S_S2_Cc_070, S_S5_Hf_107, S_NS6_Cb_042, S_S2_Cb_068, S_NS6_Am_029, S_S2_Shsn_087, S_NS1_Blss_009, S_S2_Af_119, S_S5_Bg_092_S37, S_S2_El_072, S_S2_El_073_S31, S_NS6_Pm_056 </details> 

```
Alignment sequence length ....................: `18,489 `
Version ......................................: `FastTree Version 2.1.11 Double precision (No SSE3)` 
Alignment ....................................: standard input 
Info .........................................: Amino acid distances: BLOSUM45, Joins: balanced, Support: SH-like 1000 
Search .......................................: Normal +NNI +SPR (2 rounds range 10) +ML-NNI opt-each=1 **TopHits** ......................................: 1.00\*sqrtN, close=default, refresh=0.80 
ML Model .....................................: Jones-Taylor-Thornton, CAT approximation with 20 rate categories
Info .........................................: Ignored unknown character X (seen 22,680 times) 
Refining topology ............................: 27 rounds ME-NNIs, 2 rounds ME-SPRs, 14 rounds ML-NNIs 
Branch length ................................: Total = 8.278 after 21.78 sec 
ML-NNI round 1 ...............................: LogLk = -658939.021, NNIs = 7, max delta = 100.73, Time = 65.58 sec 
ML-NNI round 2 ...............................: LogLk = -608967.117, NNIs = 3, max delta = 10.95, Time = 95.80 sec 
ML-NNI round 3 ...............................: LogLk = -608959.059, NNIs = 0, max delta = 0.00, Time = 100.13 sec 
Final ML-NNI round .........................: LogLk = -608913.064, NNIs = 0, Time = 133.33 sec
Optimize all lengths .........................: LogLk = -608912.314, Time = 143.42 sec
```
### July 15th, 2025

With my R code, I was able to determine which samples are repeated between w5 and w6, and w9 and w10. Which are the following:

`$List1_vs_List2
 [1] "036R_R_M"  "073R_0_M"  "1030R_L_M" "1067R_R_M" "1121R_0_M" "1151R_R_M" "1193R_R_F" "2404R_R_F" "2446R_L_F" "2529R_R_F"
[11] "3795R_0_F" "3909R_0_M" "3918R_R_M" "3969R_0_F" "3971R_L_F" "3975R_R_F"`

However, because the codes are a mess, I will have to manually check some of them. For instance, in some cases the 'F' is before the community consortia, and sometimes the sex is indicated as "Female" or "Male" right before the Illumina sample code. Some of them have the word "Pellet" and some others don't, so I have to remove that word too. 

It seems like all 16 samples from week 6 are repeated in week 5. Therefore, I'll remove all sample 6 from the analyses. No repetitions are present between w9 and w10


### July 16th, 2025
I decided to re-structure the original code that I had for data formatting and manipulation. That way I had a cleaner code that shows me step by step what I was doing. This new file was derived from the original `FemMicro16S_db2.R`, and it is titled `phyloseq_t1d_db2.R`. I finalized working on it, and I finally produced 4 count tables of taxa and samples divided by consortia. I have also removed the 16 repeated samples from week 6, so now I only have w5, w7, w9, w10 and the four consortia. The four otu_tables saved as .csv files are included here too, and named: ```ps_ns1_final.csv, ps_ns6_final.csv, ps_s2_final.csv, ps_s5_final.csv```

### July 17th, 2025
I have been making progress with this 2nd PICRUSt2 approach we have decided to move on with. For reference, this is what I need for my 'local files':
1. A full-length 16S (not v3v4) alignment from my reference database with a name that will be the same for the other files but with different extensions (i.e. ns1_local.fasta)
2. A phylogenetic tree from this alignment (i.e. ns1_local.tre)
3. An hmm profile from these references (i.e. ns1_local.hmm)
4. A model from the phylogenetic tree that was built with these references (i.e. ns1_local.model)

### July 18th, 2025
I was able to write a code that takes the names of the 16S headers from my fasta file, and then creates the annotated version of the KOs and ECs with those names. So, originally I was running the code `eggnog_to_picrust.py` but the name of the assemblies were exactly the same as the genomes I was running. Because PICRUSt2 needs the names to match, I had to come up with a code that I've named `adding_genomes_ko_ec.py`. This code looks at the headers from the 16S db file provided (in my case there's a total of 209 sequences from barrnap or Sanger from 116 strains) and matches the ones from my annotated EC or KO files. Once they match, it copies the entire row of annotated KOs and ECs into a new file that is compatible with PICRUSt2.

I have added the files that are produced by this script; their names are `KO_for_picrust2.tsv` and `EC_for_picrust2.tsv`. Both of them are found in my local computer in this path: `/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2_june232025/`

### July 20th, 2025
I wrote a Python script that takes the information from my FINAL master file with the final 16S count number, and then it creates a new .txt file that has the 16S database2 ID name and the 16S count that is linked to that file. So, because I originally have 116 genomes, but a total of 209 16S entries for my database, I had to link the number of 16S from the genome file to the 16S db file, having repeated CNV numbers for its correspoinding 16S id. For instance, if I have a full-length 16S Sanger sequencing from A. muciniphila, and a full-length 16S in silico genome-extracted, then I will have two entries of the same organism with the same 16S CNV that I have to retrieve from my master file. This new file will be used for PICRUSt2 and its named `16S.txt`

My initial code `16S_CNV_IMG.py` took the original or alternative names found for the original strain or its closest sister taxa from the rrnDB-v5.1 database and printed out the median, count, min, max of 16S CNV. Then I took the median as the number to use for the final 16S CNV. After that, the code that I used to do what I explained in the paragraph above is called `adding_16S_CNV.py`

### July 21st, 2025
I have moved on with the PICRUSt2 analysis. I have created in my local computer the following reference files and input files for each consortia. For simplicity, I am only including the NS1 example:

<ins>Working directory</ins>:
`/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2_june232025`

<ins>Reference files</ins> (only references! Not ASVs from the FemMicro output!):
- `ns1_local_file/ns1_local_file.fasta` this one has the 16S full-length reference sequences (v1v9 and barrnap) aligned with MAFFT.
- `ns1_local_file/ns1_local_file.hmm` this one has the hmm file created from the previous file. It needs to be done from the alignment, not from the unaligned reference sequences
- `ns1_local_file/ns1_local_file.tre` this one is the tree created with RAxML from the aligned reference file
- `ns1_local_file/ns1_local_file.model` this one is the model used to create said tree.

<ins>General annotation files that will be used for all 4 consortia</ins>:
- `picrust2_formatted_annotations/16S.txt` This one has the 16S CNV for all the genomes in all consortia that I generated with the `16S_CNV_IMG.py` script followed by `adding_16S_CNVs.py`
- `picrust2_formatted_annotations/EC_for_picrust2.tsv` This one is the ECs for all the genomes in all consortia that I generated with `eggnog_to_picrust2.py` followed by `adding_genomes_ko_ec.py`
- `picrust2_formatted_annotations/KO_for_picrust2.tsv` This one is the KOs for all the genomes in all consortia that I generated with `eggnog_to_picrust2.py` followed by `adding_genomes_ko_ec.py`

<ins>Input files</ins>: 
- `ns1_input_files/ns1_ASVs.fasta` the unaligned fasta sequences from the FemMicro output.
- `ps_ns1_asv_final_renamed.csv` This file was used to generate the biom file below. This one CANNOT be used directly as input for PICRUSt2. This file is the ASV count table from FemMicro16S. Each ASV was renamed to "ASV1","ASV2","ASV3"..."ASV85". 
- `ps_ns1_asvs.biom` The biom file formatted in an OTU table, this was generated from my R script and then formatted properly with biom (see below)

<ins>Formatting input biom file</ins>:
To convert from `ps_ns1_asv_final_renamed.csv` to biom, I first had to transform the .csv file into a .tsv file by simply doing  `sed 's/,/\t/g' ps_ns1_asv_final_renamed.csv > ps_ns1_asvs.tsv`. Then, I added the first 'cell' of the first row to be `#OTU ID`. After that, I replaced the "" that are surrounding the sample names and ASV IDs (this should be fixable in my R code...) by going into vi, and simply putting `%s/"//g`. Finally I converted my .tsv file into biom by doing `biom convert -i ps_ns1_asvs.tsv -o ps_ns1_asvs.biom --to-hdf5 --table-type="OTU table"`
- ps_ns1_asv_final_renamed.csv` the ASV count table from FemMicro. Each ASV was renamed to 'ASV1, ASV2, ASV3...'

<ins>**Placement sequence command run**</ins>: `place_seqs.py -s ns1_input_files/ns1_ASVs.fasta --ref_dir ns1_local_file/ -o ns1_output/ns1_placed_seqs.tre -p 10 --intermediate ns1_output/placement_working_ns1`

<ins>**Hidden state prediction commands**</ins>: 
`hsp.py -t ns1_output/ns1_placed_seqs.tre --observed_trait_table picrust2_formatted_annotations/16S.txt -p 10 -n -o ns1_output/ns1_16S_nsti.predicted.tsv`
`hsp.py -t ns1_output/ns1_placed_seqs.tre --observed_trait_table picrust2_formatted_annotations/EC_for_picrust2.tsv -p 10 -n -o ns1_output/ns1_EC_nsti.predicted.tsv`
`hsp.py -t ns1_output/ns1_placed_seqs.tre --observed_trait_table picrust2_formatted_annotations/KO_for_picrust2.tsv -p 10 -n -o ns1_output/ns1_KO_nsti.predicted.tsv`

<ins>**Metagenome prediction commands**</ins>:
`metagenome_pipeline.py -i ns1_input_files/ps_ns1_asvs.biom -m ns1_output/ns1_16S_nsti_predicted.tsv -f ns1_output/ns1_KO_nsti.predicted.tsv -o ns1_output/ns1_KO_metagenome_out --wide_table --strat_out`
`metagenome_pipeline.py -i ns1_input_files/ps_ns1_asvs.biom -m ns1_output/ns1_16S_nsti_predicted.tsv -f ns1_output/ns1_EC_nsti.predicted.tsv -o ns1_output/ns1_EC_metagenome_out --wide_table --strat_out`

<ins>**Pathway prediction commands**</ins>L
`pathway_pipeline.py -i ns1_output/ns1_EC_metagenome_out/pred_metagenome_strat.tsv -m pathway_mapfiles/metacyc_path2rxn_struc_filt_pro.txt -o ns1_output/ns1_pathway_out --intermediate ns1_output/ns1_EC_intermediate_pathway_files --coverage -p 10 --per_sequence_contrib --per_sequence_abun ns1_output/ns1_EC_metagenome_out/seqtab_norm.tsv --per_sequence_function ns1_output/ns1_EC_metagenome_out/pred_metagenome_strat.tsv --wide_table`







