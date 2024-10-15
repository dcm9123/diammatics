import argparse

parser = argparse.ArgumentParser(description="This program was created to generate quast jobs that can be run through a slurm file. Simply copy and paste the name of the output file into the slurm file to run it.")
print('''This job requires the following information:
        --input_path [Path to where the genome assembly directories are]
        --prefix_dir [Name of the prefix of the directories containing the input files [i.e. unicycler_test]
        --output_file [Name of the bash output file to have the quast commands written (i.e. quast_job_1.sh)
        --threads [NUM]
        --file_extension [STR] [i.e. .fasta, .fa, .gff]
        --output_dir [STR] [prefix to use for the output quast directory [i.e. quast_output]
        --sample_file [STR] [name of the file where the sample names are located [i.e. list_files.txt]
        ''')

parser.add_argument('--input_path', type=str)
parser.add_argument('--output_file', type=str)
parser.add_argument('--threads', type=int)
parser.add_argument('--file_extension', type=str)
parser.add_argument('--output_dir', type=str)
parser.add_argument('--sample_file', type=str)
parser.add_argument('--prefix_dir',type=str)

args = parser.parse_args()

input_path = args.input_path
file_out = args.output_file
threads = args.threads
file_extension = args.file_extension
output_dir = args.output_dir
sample_file = args.sample_file
prefix_dir = args.prefix_dir

output_file = open(file_out,"w")
output_file.write("#!/bin/bash\n")
with open(sample_file,'r') as file:
    lines = file.readlines()
    for line in lines:
        line = line.strip()
        output_file.write("quast --output-dir "+output_dir+"_"+line+" --threads "+str(threads)+" "+input_path+prefix_dir+"_"+line+"/assembly."+file_extension+"\n")
