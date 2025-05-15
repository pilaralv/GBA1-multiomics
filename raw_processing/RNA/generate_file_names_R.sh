# Initialize an empty variable to store the file names string
file_names_str=""

# Loop through the sample names file and append each name to the string
while IFS= read -r line; do
  file_names_str+="\"./${line}"/"${line}.genes.csv\", "  # Assuming the file extension is .csv
done < sample_names_ips_salmon.txt

# Remove the trailing comma and space at the end
file_names_str=${file_names_str%, *}

# Print the generated string
echo "file_names <- c(${file_names_str})"