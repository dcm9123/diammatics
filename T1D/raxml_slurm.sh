#!/bin/bash

#SBATCH --partition=synergy,cpu2019,cpu2021,cpu2022,cpu2023,cpu2017-bf05,cpu2019-bf05
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00
#SBATCH --mem=128G
#SBATCH --error=raxml.%J.err
#SBATCH --output=raxml.%J.out

raxml-ng --bootstrap --msa /bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/16S_data/16S_after_fixing_and_discussion/good_files_to_use/MAFFT_consortia/rm_dupt1d_db2_withDups_S_S2_species_aligned.fasta --prefix S2 --seed 100 --threads 2 --model GTR+G

raxml-ng --bootstrap --msa /bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/16S_data/16S_after_fixing_and_discussion/good_files_to_use/MAFFT_consortia/rm_dupt1d_db2_withDups_S_S5_species_aligned.fasta --prefix S5 --seed 100 --threads 2 --model GTR+G

raxml-ng --bootstrap --msa /bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/16S_data/16S_after_fixing_and_discussion/good_files_to_use/MAFFT_consortia/rm_dupt1d_db2_withDups_S_NS6_species_aligned.fasta --prefix NS6 --seed 100 --threads 2 --model GTR+G

raxml-ng --search --msa /bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/16S_data/16S_after_fixing_and_discussion/good_files_to_use/MAFFT_consortia/rm_dupt1d_db2_withDups_S_S2_species_aligned.fasta --model GTR+G --prefix S2_ML --seed 100 threads 1 --force

raxml-ng --search --msa /bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/16S_data/16S_after_fixing_and_discussion/good_files_to_use/MAFFT_consortia/rm_dupt1d_db2_withDups_S_S5_species_aligned.fasta --model GTR+G --prefix S5_ML --seed 100 threads 1 --force

raxml-ng --search --msa /bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/16S_data/16S_after_fixing_and_discussion/good_files_to_use/MAFFT_consortia/rm_dupt1d_db2_withDups_S_NS6_species_aligned.fasta --model GTR+G --prefix NS6_ML --seed 100 threads 1 --force
