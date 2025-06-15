#This script highlights the basics of bwa indexing and aligning

bwa mem <ref.fasta> -t 32 <R1.fastq> <R2.fastq> > <sam_output.sam>

#Simple conversion from bam to sam, telling samtools that the input is sam.file (-S) and the output is .bam (-b)
samtools view -b -S <file.sam> > <file.bam>

#This sorts the .bam file (arranging the alignments from the left most position)
samtools sort -@ 16 input.bam -o sorted_output.bam

#Indexing the sorted .bam file (creates a table of contents for fast access to specific alignment positions)
samtools index -@ 16 sorted.bam

#Getting average depth of coverage per sample
samtools depth output.sorted.bam | awk '{sum+=$3} END { print "Average coverage = ",sum/NR}'

#This can be easily changed for multiple files by doing this, which will simply print the number on the command line for the sample type
for file in *.sorted_indexed.bam; do samtools depth ${file} | awk -v f=$file '{sum = sum+$3} END { print f " = ", sum/NR}'; done
