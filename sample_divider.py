"""
Author: Daniel Castaneda Mogollon
Date: 08/17/2023
Purpose: This script will divide a file that contains the sample names (typically prefix of illumina files) into smaller chunks of files given a requested size from the user.
"""

import argparse
import math

parser = argparse.ArgumentParser(description="This script will divide a file that contains the sample names (typically prefix of illumina files) into smaller chuncks of files given a requested size from the user")
print('''--sample_file [STR] [The name of the file that contains all the sample names for that specific folder)
--chunks [INT] [The number of chuncks that you wish to divide each file] (i.e. 10 will make a 10 samples per file)
      ''')

parser.add_argument('--sample_file',type=str, help="The name of the file that contains all the sample names for that specific folder")
parser.add_argument('--chunks',type=int, help="The number of chuncks that you wish to divide each file (i.e. 10 will make a 10 samples per file)")
args = parser.parse_args()
sample_file = args.sample_file
chunks = args.chunks

with open(sample_file,"r") as file:                                         #Opens the sample file provided by the user
    lines = file.readlines()                                                #Stores the lines into the lines variable
    total_lines = len(lines)                                                #Counts and stores the number of lines
    piece = total_lines/chunks                                              #Divides the number of lines by the number of chunks provided
    piece = math.ceil(piece)                                                #Rounds up the size of lines in each output file
    for i in range(0,total_lines,piece):                                    #For loop from 0 to number of lines and by size of each piece
        chunk_lines = lines[i:i + piece]                                    #Stores the lines from in chunks for each file.
        with open(str(i//(piece)+1)+"_chunk_"+sample_file, "w") as chunk_file:  #Opens individual output files by a counter.
            chunk_file.writelines(chunk_lines)                              #Writes the n pieces of lines in each file
