for file in seq*.fasta; do
  # Extract number
  num=${file#seq}
  num=${num%.fasta}

  # Format number with leading zeros to 3 digits
  newnum=$(printf "%03d" "$num")

  # Create new filename
  newfile="seq${newnum}.fasta"

  # Rename file if needed
  if [[ "$file" != "$newfile" ]]; then
    mv "$file" "$newfile"
  fi
done
