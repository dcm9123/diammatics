"""
Author: Daniel Castaneda Mogollon
Date: 10/08/2024
Purpose: To print a series of bash commands that will generate fastqc output
"""

import argparse

parser = argparse.ArgumentParser(description="This program is intended to generate a series of fastqc commands that ban be run through a slurm file")
print('''This job requires the following parameters/information:
        --input_path [Path to where the polished or raw files are located]
        --output_path [Path to where the fastqc files will be deposited]
        --sample_file [A txt file that shows the name of the files without their extension and R1/R2 suffix, i.e. Sample1_S01]
        --output_file [The file name where the fastqc commands will be written, make sure to provide the .sh extension]
        ''')

parser.add_argument('--input_path', type = str)
parser.add_argument('--output_path', type = str)
parser.add_argument('--sample_file', type = str)
parser.add_argument('--output_file', type = str)

args = parser.parse_args()
input_path = args.input_path
output_path = args.output_path
sample_file = args.sample_file
output_file = args.output_file

out_file = open(output_file, "w")
out_file.write("#!/bin/bash\n")
with open(sample_file, "r") as file:
    lines = file.readlines()
    for line in lines:
        line = line.strip()
        out_file.write("fastqc -o "+output_path+" "+input_path+line+"_polished_1.fastq "+input_path+line+"_polished_2.fastq\n")
