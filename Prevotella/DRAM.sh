!#/bin/bash
#DRAM.py annotate -i unicycler_output/unicycler_output_S_32_S18/assembly.fasta -o DRAM_output/S_32_S18 --use_uniref --threads 50
DRAM.py distill -i DRAM_output/S_32_S18/annotations.tsv -o DRAM_output/S_32_S18/genome_summaries --trna_path DRAM_output/S_32_S18/trnas.tsv --rrna_path DRAM_output/S_32_S18/rrnas.tsv
#DRAM.py annotate -i unicycler_output/unicycler_output_S_FM_ATCC_S22/assembly.fasta -o DRAM_output/S_FM_ATCC_S22 --use_uniref --threads 50
DRAM.py distill -i DRAM_output/S_FM_ATCC_S22/annotations.tsv -o DRAM_output/S_FM_ATCC_S22/genome_summaries --trna_path DRAM_output/S_FM_ATCC_S22/trnas.tsv --rrna_path DRAM_output/S_FM_ATCC_S22/rrnas.tsv
#DRAM.py annotate -i unicycler_output/unicycler_output_S_FM_DSM_S23/assembly.fasta -o DRAM_output/S_FM_DSM_S23 --use_uniref --threads 50
DRAM.py distill -i DRAM_output/S_FM_DSM_S23/annotations.tsv -o DRAM_output/S_FM_DSM_S23/genome_summaries --trna_path DRAM_output/S_FM_DSM_S23/trnas.tsv --rrna_path DRAM_output/S_FM_DSM_S23/rrnas.tsv
#DRAM.py annotate -i unicycler_output/unicycler_output_S_IMP_BC_01_S43/assembly.fasta -o DRAM_output/S_IMP_BC_01_S43 --use_uniref --threads 50
DRAM.py distill -i DRAM_output/S_IMP_BC_01_S43/annotations.tsv -o DRAM_output/S_IMP_BC_01_S43/genome_summaries --trna_path DRAM_output/S_IMP_BC_01_S43/trnas.tsv --rrna_path DRAM_output/S_IMP_BC_01_S43/rrnas.tsv
#DRAM.py annotate -i unicycler_output/unicycler_output_S_IMP_BC_02_S44/assembly.fasta -o DRAM_output/S_IMP_BC_02_S44 --use_uniref --threads 50
DRAM.py distill -i DRAM_output/S_IMP_BC_02_S44/annotations.tsv -o DRAM_output/S_IMP_BC_02_S44/genome_summaries --trna_path DRAM_output/S_IMP_BC_02_S44/trnas.tsv --rrna_path DRAM_output/S_IMP_BC_02_S44/rrnas.tsv
#DRAM.py annotate -i unicycler_output/unicycler_output_S_IMP_BC_04_S45/assembly.fasta -o DRAM_output/S_IMP_BC_04_S45 --use_uniref --threads 50
DRAM.py distill -i DRAM_output/S_IMP_BC_04_S45/annotations.tsv -o DRAM_output/S_IMP_BC_04_S45/genome_summaries --trna_path DRAM_output/S_IMP_BC_04_S45/trnas.tsv --rrna_path DRAM_output/S_IMP_BC_04_S45/rrnas.tsv
#DRAM.py annotate -i unicycler_output/unicycler_output_S_IMP_PBS-GLY_S42/assembly.fasta -o DRAM_output/S_IMP_PBS-GLY_S42 --use_uniref --threads 50
DRAM.py distill -i DRAM_output/S_IMP_PBS-GLY_S42/annotations.tsv -o DRAM_output/S_IMP_PBS-GLY_S42/genome_summaries --trna_path DRAM_output/S_IMP_PBS-GLY_S42/trnas.tsv --rrna_path DRAM_output/S_IMP_PBS-GLY_S42/rrnas.tsv
#DRAM.py annotate -i unicycler_output/unicycler_output_S_Josh_DNF00844_S49/assembly.fasta -o DRAM_output/S_Josh_DNF00844_S49 --use_uniref --threads 50
DRAM.py distill -i DRAM_output/S_Josh_DNF00844_S49/annotations.tsv -o DRAM_output/S_Josh_DNF00844_S49/genome_summaries --trna_path DRAM_output/S_Josh_DNF00844_S49/trnas.tsv --rrna_path DRAM_output/S_Josh_DNF00844_S49/rrnas.tsv
#DRAM.py annotate -i unicycler_output/unicycler_output_S_Josh_KA00225_S25/assembly.fasta -o DRAM_output/S_Josh_KA00225_S25 --use_uniref --threads 50
DRAM.py distill -i DRAM_output/S_Josh_KA00225_S25/annotations.tsv -o DRAM_output/S_Josh_KA00225_S25/genome_summaries --trna_path DRAM_output/S_Josh_KA00225_S25/trnas.tsv --rrna_path DRAM_output/S_Josh_KA00225_S25/rrnas.tsv
#DRAM.py annotate -i unicycler_output/unicycler_output_S_Josh_KA00625_S24/assembly.fasta -o DRAM_output/S_Josh_KA00625_S24 --use_uniref --threads 50
DRAM.py distill -i DRAM_output/S_Josh_KA00625_S24/annotations.tsv -o DRAM_output/S_Josh_KA00625_S24/genome_summaries --trna_path DRAM_output/S_Josh_KA00625_S24/trnas.tsv --rrna_path DRAM_output/S_Josh_KA00625_S24/rrnas.tsv
