#!/usr/bin/env bash


# INPUT
file_lit="${snakemake_input[0]}"
#study=INTERVAL
study=CHRIS
study_id_col_name="study_id"
chr_col_name="chr"
pos_col_name="POS"
snpid_col_name="SNPID"  # Adding SNPID column name

gwas_dir=/exchange/healthds/pQTL/results/${study}/qced_sumstats/outputs/

# OUTPUT
out_dir=$TMPDIR
base_name=$(basename "$file_lit" .csv)
file_lit_SEQID_SNPID=${out_dir}/file_lit_SEQID_SNPID_${study}
out_file="${snakemake_output[0]}"
out_nofile=${out_dir}/out_nofile_${study}.tsv
# Clean option
clean=true  # Set this to true to enable cleanup


# START

# Extract and display the header
header=$(head -n 1 "$file_lit")
echo "Header: $header"

# Function to get the column index for a given column name
get_col_index() {
    col_name=$1
    echo "$header" | awk -F, -v col_name="$col_name" '{
        for(i=1; i<=NF; i++) {
            gsub(/\"/, "", $i); # Remove quotes
            if($i==col_name) {
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
snpid_col=$(get_col_index "$snpid_col_name")  # Get SNPID column index

# Verify if the column indices are correctly identified
echo "$study_id_col_name column index: $study_id_col"
echo "$chr_col_name column index: $chr_col"
echo "$pos_col_name column index: $pos_col"
echo "$snpid_col_name column index: $snpid_col"

if [ ! -f $file_lit_SEQID_SNPID ]; then
        echo "Creating initial File!"
        # extract the columns: new file needs seqid in first column, snpid in second column, and chr:pos1-pos2 in third column
        awk -F, -v OFS='\t' -v study_id_col="$study_id_col" -v chr_col="$chr_col" -v pos_col="$pos_col" -v snpid_col="$snpid_col" -v gwas_dir="$gwas_dir" '
        NR > 1 {
            chr=$chr_col;
            pos=$pos_col;
            snpid=$snpid_col;
            chr_pos=chr ":" pos "-" pos;
            gsub(/\"/, "", $study_id_col); # Remove quotes from study_id
            gsub(/\"/, "", chr_pos);       # Remove quotes from chr_pos
            gsub(/\"/, "", snpid);         # Remove quotes from snpid
            gwas_path=gwas_dir $study_id_col "/" $study_id_col ".gwaslab.tsv.gz";
            print $study_id_col, snpid, chr_pos, gwas_path;
        }' "$file_lit" > "$file_lit_SEQID_SNPID"
fi

# Append the header to out_file with additional columns "SNPID_lit" and "filename"
first_gwas_path=$(awk 'NR==1{print $4}' "$file_lit_SEQID_SNPID")
first_gwas_header=$(zcat -f "$first_gwas_path" | head -n 1)

# Append the header to out_file with additional columns "SNPID_lit", "filename", and "study"
echo -e "SNPID_lit\tfilename\tstudy\t$first_gwas_header" > "$out_file"
echo -e "SNPID_lit\tfilename\tstudy" > "$out_nofile"

cat $file_lit_SEQID_SNPID | { cat ; echo ; } | while read line; do
    seqid=$(echo $line | awk '{print $1}')
    snpid=$(echo $line | awk '{print $2}')
    chr_pos_pos=$(echo $line | awk '{print $3}')
    filename=$(echo $line | awk '{print $4}')
    if [ -e $filename ]; then
        echo "file exists"
        tabix_output=$(tabix $filename $chr_pos_pos)
        matched=false
        echo "$tabix_output" | while read tabixout; do
            if [ -n "$tabixout" ]; then
                # Check if any field in the tabix output matches the snpid
                if echo "$tabixout" | awk -v snpid="$snpid" '{for(i=1;i<=NF;i++) if($i==snpid) exit 0; exit 1}'; then
                    paste <(echo $snpid) <(echo "$filename") <(echo "$study") <(echo "$tabixout") >> $out_file
                    matched=true
                    break
                fi
            fi
        done
        if [ "$matched" = false ]; then
            echo "No matching SNPID found for $snpid in $filename"
        fi
    else
        echo "file $filename does not exist"
        echo -e "$seqid\t$filename\t$study" >> $out_nofile
    fi
done

# CHECKS:
# Check the row counts
input_row_count=$(($(wc -l < "$file_lit") - 1))  # Exclude header
output_row_count=$(($(wc -l < "$out_file") - 1))  # Exclude header

if [ "$input_row_count" -ne "$output_row_count" ]; then
    echo "Warning: The number of rows in $file_lit ($input_row_count) does not match the number of rows in $out_file ($output_row_count)"
else
    echo "The number of rows match between $file_lit and $out_file"
    if [ "$clean" = true ]; then
        echo "Cleaning up intermediate files..."
        rm -f "$file_lit_SEQID_SNPID" "$out_nofile"
    fi
fi
