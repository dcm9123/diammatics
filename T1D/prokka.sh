!/bin/bash

for dir in *; do prokka --outdir prokka_${dir} --prefix ${dir} --cpus 56 --kingdom Bacteria --rfam ${dir}/S*.fasta; done;
