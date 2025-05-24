from Bio import SeqIO
import os

path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/t1d_db_fixed_discussed/good_files_to_use/"
os.chdir(path)
file1 = "t1d_db_fullTaxo.fasta"
file2 = "t1d_db_species.fasta"

def reading_files(fasta_file, sample_id, fasta_preffix):
    fasta_out = open("t1d_db_updated_"+sample_id+fasta_preffix, "w")
    in_file = SeqIO.parse(fasta_file, "fasta")
    i = 0
    for record in in_file:
        if sample_id in record.id:
            i = i + 1
            SeqIO.write(record, fasta_out, "fasta")
            #fasta_out.write(re)
    print(f"There are a total of {i} records written to {sample_id}{fasta_preffix}")

def inspecting(fasta_file):
    in_file = SeqIO.parse(fasta_file, "fasta")
    for record in in_file:
        if "-" in record.seq:
            print(fasta_file)
            print(record.id)

def main():
    reading_files(file1, "S_NS1_", "fullTaxo.fasta")
    reading_files(file1, "S_NS6_", "fullTaxo.fasta")
    reading_files(file1, "S_S2_", "fullTaxo.fasta")
    reading_files(file1, "S_S5_", "fullTaxo.fasta")
    reading_files(file2, "S_NS1_", "species.fasta")
    reading_files(file2, "S_NS6_", "species.fasta")
    reading_files(file2, "S_S2_", "species.fasta")
    reading_files(file2, "S_S5_", "species.fasta")
    #inspecting(file1)

main()
