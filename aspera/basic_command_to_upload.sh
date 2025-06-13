#Daniel Castaneda Mogollon, PhD
#June 12th, 2025
#This simple script is to upload .fasta files to my username in ncbi for later submission

ascp -i Desktop/Keys/aspera.openssh -QT -l100m -d Desktop/Sycuro/Projects/Diabetes/Genomes_June2025/all_genomes_1st_run/ subasp@upload.ncbi.nlm.nih.gov:uploads/daniel.castanedamogo_ucalgary.ca_UoWD2Ptv|

#The -i is wher ethe key is
#The -Q is to leave it in quiet mode
#THE -T disables encryption
#The -l100m limits the upload to 100 megabits/sec (~12MB/sec)
#The -d is the folder with all of my genomes
