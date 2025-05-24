for file in seq*.fasta; do
  # Extract just the number part using parameter expansion and pattern matching
  num=${file#seq}
  num=${num%.fasta}

  # Count the number of digits
  case ${#num} in
    1)
      mv "$file" one_digit/
      ;;
    2)
      mv "$file" two_digits/
      ;;
    3)
      mv "$file" three_digits/
      ;;
    *)
      echo "Skipping $file (unexpected number format)"
      ;;
  esac
done
