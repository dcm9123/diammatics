"""
Author: Daniel Castaneda Mogollon
Date: 08/01/2023
Purpose: To print a list of commands that could assemble polished-reads into unicycler through a slurm script
"""


import os
import sys
import argparse

parser = argparse.ArgumentParser(description="This program was created to generate unicycler jobs that can be run through a slurm file. Simply copy and paste the output of this command into a slurm file to run its output.")
print('''This job requires the following information:
        --input_path [Path to where the sequenced files are located]
        --output_file [Name of the file for the bash job to be written]
        --threads [NUM]
        --mode [conservative,normal,bold]
        --depth_filter [FLOAT]
        --min_contig_length [FLOAT]
        --illumina_suffix [STR] [i.e. _filtered, _raw]
        --file_extension [STR] [i.e. .fastq, .fq]
        --output_dir [STR] [prefix of the output directory to where the unicycler output will be saved]
        --sample_file [STR] [name of the file where the sample names are located]
        ''')
parser.add_argument('--input_path','-i', type=str, help='Input file where the prefix of the sample names are (without the illumina sequencing suffix info)')
parser.add_argument('--output_file','-o', type=str, help='Output file where the list of commands will be generated to be included in a slurm script')
parser.add_argument('--threads','-t', type=int, help='The number of threads to be assigned to the unicycler jobs')
parser.add_argument('--mode', type=str, help='The stringency of the assembler when joining contigs (default = normal)')
parser.add_argument('--depth_filter', type=float, help='This represents the fraction of contigs lower of the chromosomal depth to be filtered out (default = 0.25)')
parser.add_argument('--min_contig_length', type=int, help='This is the minimum size of contigs that will be included in the analysis')
parser.add_argument('--illumina_suffix', type=str, help='This is the name of the illumina info after the sample preffix, i.e. _L1_R1_001, _L1_R2_001')
parser.add_argument('--file_extension', type=str, help='This is the file extension for the sequencing files, i.e. .fastq, .fq, fastq.gz, fq.gz')
parser.add_argument('--output_dir', type=str, help='This is the output directory where Unicycler will produce the assembly output')
parser.add_argument('--sample_file',type=str, help='This is the file that has the sample names without the illumina lane information, i.e. S1_Pd')
args = parser.parse_args()

num_threads = args.threads
input_path = args.input_path
output_file = args.output_file
mode = args.mode
depth_filter = args.depth_filter
min_contig_length = args.min_contig_length
illumina_suffix = args.illumina_suffix
file_extension = args.file_extension
sample_file = args.sample_file
file_out = args.output_file
dir_out = args.output_dir

output_file = open(file_out,"w")
output_file.write("#!/bin/bash\n")
with open(sample_file,'r') as file:
    lines = file.readlines()
    for line in lines:
        line = line.strip()
        output_file.write("unicycler -1 "+input_path+str(line)+illumina_suffix+"_1"+file_extension+" -2 "+input_path+str(line)+illumina_suffix+"_2"+file_extension+" --mode "+mode+" --threads "+str(num_threads)+" --keep 2 --depth_filter "+str(depth_filter)+" --min_fasta_length "+str(min_contig_length)+" --out "+dir_out+"_"+line+"\n")
