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
6. Run the contigs workflow like this: `anvi-run-workflow -w contigs -c config_default.json`. In my case, with 116 samples, I had to run the interactive node 3 times for 5 hours each. However, when rerunning it, anvi'o is a bit picky on rewriting some of the files and pick up where it left off. So, in that case I had to run this command `anvi-run-workflow - contigs -c config_default.json --additional-params --rerun-incomplete --keep-going`. This gave me the contigs database I need for downstream steps.

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

```
Input aligment file path .....................: /bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/anvio/all_genomes_NCBI/final_concatenated_proteins.fa    
Output file path .............................: /bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/anvio/all_genomes_NCBI/NCBI_phylogenomics_tree.txt       
Alignment names ..............................: S_S5_Rf_112, S_S5_Cc_098, S_S5_Sc_120, S_NS6_Hh_054, S_NS1_Fp_020_S27, S_S5_Bv_093, S_S5_Cb_096, S_S5_Pd_110, S_S5_Cc_097, S_S2_Blss_067_S30, S_S2_Cc_069, S_S5_Ef_104, S_NS6_Fp_052, S_NS6_El_044, S_NS6_Ef_050,
S_S2_Fp_079, S_NS1_Bs_004, S_NS1_Bv_007, S_NS1_Ss_026, S_NS1_Fp_019, S_S5_Rg_113, S_NS1_Bu_006, S_NS1_Ao_003,
S_S2_Bs_063, S_S2_Bv_065, S_NS1_Ca_028, S_NS1_Shsn_025, S_NS1_Bp_010, S_NS1_Em_018, S_NS6_Ef_049, S_S5_Em_101,
S_S2_Et_074_S32, S_NS1_Pd_022, S_S5_Dg_099, S_NS6_Em_046, S_S5_Et_102, S_NS1_El_014, S_S2_Ef_078_S34,
S_NS6_Bs_033, S_NS1_Sw_027, S_NS6_Bs_032, S_NS6_Blss_040, S_NS1_Ig_115, S_NS1_Am_001, S_NS6_El_045,
S_NS6_Ao_031, S_S2_Ea_077, S_S5_Am_089, S_S5_Pd_111, S_NS1_Cb_011, S_NS1_Et_015, S_NS1_He_021, S_S5_El_100,
S_NS6_Bu_117, S_S5_Ea_103, S_S2_He_082_S35, S_S2_Bd_066, S_S2_Ca_088_S36, S_S5_Af_090, S_NS1_Eh_016,
S_S2_Rf_084, S_S2_Rg_085, S_NS1_Af_002, S_NS6_Pp_057_S29, S_NS1_Pm_023, S_NS1_Bl_116, S_NS6_Bu_036,
S_S2_Am_061, S_S2_Et_075_S33, S_NS6_Bx_038_S28, S_S5_Af_091, S_NS6_He_053, S_NS6_Sw_060, S_NS6_Et_048,
S_NS6_Ss_059, S_S2_Bd_064, S_S2_Cs_071, S_S2_Rf_118, S_S2_Fp_080, S_NS6_Cb_041, S_NS1_Ba_008, S_NS6_Bu_034,
S_NS1_Cc_012, S_S2_Pd_083, S_S5_Blss_094, S_NS1_Cp_013, S_S5_He_108, S_NS6_Bv_037, S_S5_Lp_109_S38,
S_S5_Bp_095, S_NS1_Em_017, S_S2_Ea_076, S_S5_Fp_105, S_NS6_Ba_039, S_NS6_Af_030, S_NS6_Sh_058, S_NS6_Ef_051,
S_S5_Fp_106, S_S2_Hf_081, S_NS6_Em_047, S_NS6_Pd_055, S_NS6_Cp_043, S_NS1_Bt_005, S_S2_Sc_086, S_S2_Cc_070,
S_S5_Hf_107, S_NS6_Cb_042, S_S2_Cb_068, S_NS6_Am_029, S_S2_Shsn_087, S_NS1_Blss_009, S_S2_Af_119,
S_S5_Bg_092_S37, S_S2_El_072, S_S2_El_073_S31, S_NS6_Pm_056

Alignment sequence length ....................: 18,489                                                                Version ......................................: FastTree Version 2.1.11 Double precision (No SSE3)                    Alignment ....................................: standard input                                                        Info .........................................: Amino acid distances: BLOSUM45 Joins: balanced Support: SH-like 1000  Search .......................................: Normal +NNI +SPR (2 rounds range 10) +ML-NNI opt-each=1               TopHits ......................................: 1.00*sqrtN close=default refresh=0.80                                 ML Model .....................................: Jones-Taylor-Thorton, CAT approximation with 20 rate categories       Info .........................................: Ignored unknown character X (seen 22680 times)                        Refining topology ............................: 27 rounds ME-NNIs, 2 rounds ME-SPRs, 14 rounds ML-NNIs                Info .........................................: Total branch-length 8.278 after 21.78 sec                             ML-NNI round 1 ...............................: LogLk = -658939.021 NNIs 7 max delta 100.73 Time 65.58                Info .........................................: Switched to using 20 rate categories (CAT approximation)              Info .........................................: Rate categories were divided by 0.960 so that average rate = 1.0      Info .........................................: CAT-based log-likelihoods may not be comparable across runs           Info .........................................: Use -gamma for approximate but comparable Gamma(20) log-likelihoods   ML-NNI round 2 ...............................: LogLk = -608967.117 NNIs 3 max delta 10.95 Time 95.80                 ML-NNI round 3 ...............................: LogLk = -608959.059 NNIs 0 max delta 0.00 Time 100.13                 Info .........................................: Turning off heuristics for final round of ML NNIs (converged)         ML-NNI round 4 ...............................: LogLk = -608913.064 NNIs 0 max delta 0.00 Time 133.33 (final)         Optimize all lengths .........................: LogLk = -608912.314 Time 143.42
```






