#!/usr/bin/env bash

# INPUT
file_lit="${snakemake_input[0]}"
study="${snakemake_output[0]//[^a-zA-Z0-9]/_}"
study=(${study//_/ })
study="${study[8]}"
#study=INTERVAL  # Define the study
# study=CHRIS

study_id_col_name="phenotype_id"
chr_col_name="chr"
pos_col_name="POS"
snpid_col_name="SNPID"

gwas_dir=/exchange/healthds/pQTL/results/"${study}"/qced_sumstats/outputs/

# OUTPUT
out_dir=$TMPDIR
out_file="${snakemake_output[0]}"


# START

# Extract and display the header
header=$(head -n 1 "$file_lit")  # Read the original header
echo "Header: $header" 2>> "${snakemake_log[0]}"

# Function to get the column index for a given column name
get_col_index() {
    col_name=$1
    echo "$header" | awk -v col_name="$col_name" '
    BEGIN {FS=";"}  # Set field separator to ";"
    {
        for(i=1; i<=NF; i++) {
            if($i == col_name) {
                print i;
                break;
            }
        }
    }'
}

# Get the column indices for the desired columns
study_id_col=$(get_col_index "$study_id_col_name")
chr_col=$(get_col_index "$chr_col_name")
pos_col=$(get_col_index "$pos_col_name")
snpid_col=$(get_col_index "$snpid_col_name")

# Verify if the column indices are correctly identified
echo "$study_id_col_name column index: $study_id_col"
echo "$chr_col_name column index: $chr_col"
echo "$pos_col_name column index: $pos_col"
echo "$snpid_col_name column index: $snpid_col"

# Prepare the output file by appending column headers with the study suffix
first_gwas_path="$gwas_dir/$(head -n 2 "$file_lit" | tail -n 1 | cut -d ';' -f"$study_id_col")/$(head -n 2 "$file_lit" | tail -n 1 | cut -d ';' -f"$study_id_col").gwaslab.tsv.gz"
first_gwas_header=$(zcat -f "$first_gwas_path" | head -n 1 | awk -v study="$study" '{for (i=1; i<=NF; i++) printf "%s_%s;", $i, study}')

# Write the new header with the appended columns
echo "$header;$first_gwas_header" | sed 's/;$//' > "$out_file"

# Process the input file and append data from the corresponding GWAS file
tail -n +2 "$file_lit" | while IFS=";" read -r line; do
    # Extract necessary fields
    seqid=$(echo "$line" | cut -d ';' -f"$study_id_col")
    chr=$(echo "$line" | cut -d ';' -f"$chr_col")
    pos=$(echo "$line" | cut -d ';' -f"$pos_col")
    snpid=$(echo "$line" | cut -d ';' -f"$snpid_col")

    # Set filename for the corresponding GWAS file
    filename="${gwas_dir}${seqid}/${seqid}.gwaslab.tsv.gz"

    if [ -e "$filename" ]; then
        echo "file exists: $filename"

        # Run tabix to get the row that matches both chr:pos and snpid
        tabix_output=$(tabix "$filename" "$chr:$pos-$pos" | awk -v snpid="$snpid" '$0 ~ snpid {print $0}')

        if [ -n "$tabix_output" ]; then
            # Convert tabix output to CSV format (using ";" as separator) and append to the original row
            tabix_output_csv=$(echo "$tabix_output" | tr '\t' ';')
            echo "$line;$tabix_output_csv" | sed 's/;$//' >> "$out_file"
        else
            echo "ERROR: No matching SNPID found for $snpid in $filename. Exiting."
            exit 1  # Exit the script immediately
        fi
    else
        echo "ERROR: File $filename does not exist. Exiting."
        exit 1  # Exit the script immediately
    fi
done
