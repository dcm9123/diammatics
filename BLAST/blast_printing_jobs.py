#Program created by Daniel Castaneda Mogollon, PhD
#June 9th, 2024

import argparse
from Bio import SeqIO
import os

parser = argparse.ArgumentParser(description="This program prints the commands needed to run a shell script in a slurm-structured server, requiring the user to provide a database path, a path where individual fasta files are located, and an output path")
parser.add_argument("--input_path", type=str, help="This is the path where your fasta files are located")
parser.add_argument("--database", type=str, help="This is the path to the personalized or NCBI database")
parser.add_argument("--out_file", type=str, help="This is the name of the output file where the sh commands will be printed")
parser.add_argument("--out_path", type=str, help="This is the name of your output directory for the BLAST results to be printed in")

args = parser.parse_args()
input_path = args.input_path
database = args.database
out_file = args.out_file
output_path = args.out_path

def printing_commands(input_path,database,out_file,output_path):
    with open(out_file, "a") as outfile:
        outfile.write("#!/bin/bash\n")
        for filename in os.listdir(input_path):
            if filename.endswith(".fasta"):
                query_path = os.path.join(input_path,filename)
                output_path = os.path.join(input_path,filename+".txt")
                outfile.write(f"blastn -query {query_path} -db {database} -out {output_path} -outfmt '7 salltitles score evalue nident pident qlen qstart qend sstart ssend'\n")
