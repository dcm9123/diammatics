#!/bin/bash

#SBATCH --partition=synergy,cpu2019,cpu2021,cpu2022,cpu2023,cpu2017-bf05,cpu2019-bf05
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=2-00:00:00
#SBATCH --mem=120G
#SBATCH --error=eggnog_hmmer.%J.err
#SBATCH --output=eggnog_hmmer.%J.out

for dir in /bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/genomes_unicycler/unicycler/1st_run/prokka_output/prokka_S_*; do cd $dir; emapper.py -i *.faa -d Bacteria --database /bulk/IMCshared_bulk/shared/dbs/eggnog_db/hmmer/ -o $dir --cpu 0 --itype proteins --output_dir .; cd ../; done;
