from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description = "This program will rename the fasta reacords of your file, followed by a filtering provided by the user")
parser.add_argument("--in_file", type=str, help="The name of your fasta file")
parser.add_argument("--filter_size", type=int, help="This is the size in bp of the fasta sequences you wish to keep (i.e. 1200 will keep anything 1200 and higher)")
args = parser.parse_args()
in_file = args.in_file
filter_size = args.filter_size

i=1
out_file = str(in_file)+"_filtered.fasta"
for record in SeqIO.parse(in_file,"fasta"):
    with open(out_file,"a") as output_handle:
        if(len(record.seq)>=filter_size):
            record.id=str(in_file)+"_"+str(i)
            record.description=""
            SeqIO.write(record,output_handle,"fasta")
            i=i+1

    print(f"{record.id}\t{len(record.seq)}")
