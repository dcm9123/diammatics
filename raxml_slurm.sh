#!/bin/bash

#SBATCH --partition=synergy,cpu2019,cpu2021,cpu2022,cpu2023,cpu2017-bf05,cpu2019-bf05
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00
#SBATCH --mem=64G
#SBATCH --error=raxml.%J.err
#SBATCH --output=raxml.%J.out

#raxml-ng --bootstrap --msa mafft/NS1_aligned_species.fasta --prefix NS1 --seed 100 --threads 4 --model GTR+G
raxml-ng --bootstrap --msa mafft/NS6_aligned_species.fasta --prefix NS6 --seed 100 --threads 4 --model GTR+G
#raxml-ng --bootstrap --msa mafft/S2_aligned_species.fasta --prefix S2 --seed 100 --threads 4 --model GTR+G
#raxml-ng --bootstrap --msa mafft/S5_aligned_species.fasta --prefix S5 --seed 100 --threads 4 --model GTR+G
raxml-ng --search --msa mafft/NS1_aligned_species.fasta --model GTR+G --prefix NS1_ML --seed 100 threads 2
raxml-ng --search --msa mafft/NS6_aligned_species.fasta --model GTR+G --prefix NS6_ML --seed 100 threads 2
raxml-ng --search --msa mafft/S2_aligned_species.fasta --model GTR+G --prefix S2_ML --seed 100 threads 2
raxml-ng --search --msa mafft/S5_aligned_species.fasta --model GTR+G --prefix S5_ML --seed 100 threads 2
