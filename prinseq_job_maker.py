"""
Author: Daniel Castaneda Mogollon
Date: 10/07/2024
Purpose: To print a list of commands that could polish raw reads through a slurm script
"""

import os
import argparse

parser = argparse.ArgumentParser(description="This program was created to generate a prinseq job that can be run through a slurm file.")
print('''This job requires the following parameters/information:
        --input_path [Path to where your raw .fastq files are]
        --output_path [Path to where the polished reads should be deposited]
        --output_file [Name of the file where the commands will be written]
        --threads [number of threads to run in the slurm file]
        --sample_list [The path and name of the file that contains the names of your samples without the 'R1_L001_001', just the preffix]
        --trim_left [The value to assign to prinseq to trim left]
        --trim_right [The value to assign to prinseq to trim from the right]
        --lc_method [The LC method to run in prinseq, i.e. dust] 
        --lc_threshold [i.e. 7]
        --trim_qual_type [i.e mean]
        --trim_qual_window [i.e. 10]
        --trim_qual_step [i.e. 2]
        --trim_qual_rule [i.e. lt]
        --trim_qual_left [i.e. 30]
        --trim_qual_right [i.e. 30]
        --min_length [i.e. 60]
        --nx_max_n [i.e. 15]
        ''')

parser.add_argument('--input_path', type=str)
parser.add_argument('--output_path', type=str)
parser.add_argument('--output_file', type=str)
parser.add_argument('--threads', type=int)
parser.add_argument('--sample_list', type=str)
parser.add_argument('--trim_left', type=int)
parser.add_argument('--trim_right', type=int)
parser.add_argument('--lc_method', type=str)
parser.add_argument('--lc_threshold', type=int)
parser.add_argument('--trim_qual_type', type=str)
parser.add_argument('--trim_qual_window', type=int)
parser.add_argument('--trim_qual_step', type=int)
parser.add_argument('--trim_qual_rule', type=str)
parser.add_argument('--trim_qual_left', type=int)
parser.add_argument('--trim_qual_right', type=int)
parser.add_argument('--min_len', type=int)
parser.add_argument('--ns_max_n', type=int)

args = parser.parse_args()
input_path = args.input_path
output_path = args.output_path
threads = args.threads
sample_list = args.sample_list
trim_left = args.trim_left
trim_right = args.trim_right
lc_method = args.lc_method
lc_threshold = args.lc_threshold
trim_qual_type = args.trim_qual_type
trim_qual_window = args.trim_qual_window
trim_qual_step = args.trim_qual_step
trim_qual_rule = args.trim_qual_rule
trim_qual_left = args.trim_qual_left
trim_qual_right = args.trim_qual_right
min_len = args.min_len
ns_max_n = args.ns_max_n
output_file = args.output_file

out = open(output_file, "w")
out.write("#!/bin/bash\n")
with open(sample_list,"r") as file:
    lines = file.readlines()
    for line in lines:
        line = line.strip()
        out.write("perl utils/scripts/prinseq-lite.pl -fastq "+input_path+line+"_R1_001.fastq -fastq2 "+input_path+line+"_R2_001.fastq -trim_left "+str(trim_left)+" -trim_right "+str(trim_right)+" -out_good "+output_path+line+"_polished -out_bad null -lc_method "+lc_method+" -lc_threshold "+str(lc_threshold)+" -derep 1 -trim_qual_type "+trim_qual_type+" -trim_qual_window "+str(trim_qual_window)+ " -trim_qual_step "+str(trim_qual_step)+" -trim_qual_rule "+trim_qual_rule+" -trim_qual_left "+str(trim_qual_left)+" -trim_qual_right "+str(trim_qual_right)+" -min_len "+str(min_len)+" -ns_max_n "+str(ns_max_n)+"\n")
