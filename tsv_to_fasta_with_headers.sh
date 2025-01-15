#Daniel Castaneda Mogollon, PhD & Kevin Muirhead, Msc
#January 15th, 2025
#This code was created in order to generate a fasta file from a .tsv file that has an ASV id column and its corresponding sequence

IFS=$'\n'; for row in $(tail -n+2 dada2_all_databases_merged.tsv); do echo $row; asv_id=$(echo $row | cut -d ',' -f1 | sed 's/"//g'); echo $asv_id; seq=$(echo $row | cut -d ',' -f2 | sed 's/"//g'); echo ">${asv_id}\n${seq}"; done
