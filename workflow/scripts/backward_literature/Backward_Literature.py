import pandas as pd
import sqlite3
import json
import click


def load_config(config_path):
    with open(config_path, "r") as file:
        config = json.load(file)
    return config


@click.command()
@click.option("-c", "--config_path", required=True, help="Configuration file path")
@click.option("-f", "--file_path", required=True, help="Table file path")
@click.option(
    "-t",
    "--input_type",
    required=True,
    type=click.Choice(["regional_table", "unique_signal_table"], case_sensitive=False),
    help="type of table to process",
)
@click.option("--suffix", required=True, help="Output suffix", default="_bl_ann")
@click.option("--sep", required=True, help="column separator suffix", default=";")
def backward_literature_review_annotation(config_path, file_path, input_type, suffix, sep):
    """
    Main function to orchestrate the literature review annotation process.
    Handles the connection to the database and processes the tables (regional results - i.e. locus breaker - or unique signals)
    to find matches based on SeqID and UniProt.

    Parameters:
    - config_path: str
        Path to the JSON configuration file which contains necessary settings like database path,
        reference genome version, and results table details (file paths, column names, etc.).
    - input_type: str
        Specifies the type of table to process. Can be either 'regional_table' (for results that have significant 'regions' to check against literature) or 'unique_signal_table',
        which determines the dataset being processed.

    Outputs:
    - The function updates the regional or unique signal tables with matching results and saves the updated
      CSV file.
    """
    # Load the configuration file
    config = load_config(config_path)
    # overwrite the file_path
    config[input_type]["file_path"] = file_path

    # Connect to the SQLite database
    conn = sqlite3.connect(config["database_path"])

    # Determine the position column based on the reference genome version
    reference_genome = config["reference_genome"]
    pos_column = "pos37" if reference_genome == "37" else "pos38"

    # Process data based on the input type: regional table or unique signal table
    if input_type == "regional_table":
        regional_table = pd.read_csv(config["regional_table"]["file_path"], sep=sep)
        regional_table["soma_match"] = "NO"
        regional_table["soma_matching_study"] = ""
        regional_table["soma_matching_number_ids"] = ""
        regional_table["soma_matching_signals"] = ""  # Add matching signal for SeqID
        regional_table["uniprot_match"] = "NO"
        regional_table["unip_matching_study"] = ""
        regional_table["unip_matching_number_ids"] = ""
        regional_table["unip_matching_signals"] = ""  # Add matching signal for UniProt

        start_col = config["regional_table"]["start_column"]
        end_col = config["regional_table"]["end_column"]
        protein_col = config["regional_table"]["protein_column"]
        uniprot_col = config["regional_table"]["uniprot_column"]  # Use uniprot_column

        process_table(
            conn,
            regional_table,
            start_col,
            end_col,
            protein_col,
            uniprot_col,
            pos_column,
        )

        output_path = config["regional_table"]["file_path"].replace(".csv", f"{suffix}.csv")
        regional_table.to_csv(output_path, index=False, sep=sep)

    elif input_type == "unique_signal_table":
        unique_signal_table = pd.read_csv(config["unique_signal_table"]["file_path"])
        unique_signal_table["soma_match"] = "NO"
        unique_signal_table["soma_matching_study"] = ""
        unique_signal_table["soma_matching_number_ids"] = ""
        unique_signal_table["soma_matching_signals"] = ""  # Add matching signal for SeqID
        unique_signal_table["uniprot_match"] = "NO"
        unique_signal_table["unip_matching_study"] = ""
        unique_signal_table["unip_matching_number_ids"] = ""
        unique_signal_table["unip_matching_signals"] = ""  # Add matching signal for UniProt

        chr_col = config["unique_signal_table"]["chr_column"]
        pos_col = config["unique_signal_table"]["pos_column"]
        protein_col = config["unique_signal_table"]["protein_column"]
        uniprot_col = config["unique_signal_table"]["uniprot_column"]  # Use uniprot_column
        region_size = config["unique_signal_table"]["region_size"]

        for index, row in unique_signal_table.iterrows():
            start = row[pos_col] - region_size
            end = row[pos_col] + region_size
            process_table_row(
                conn,
                unique_signal_table,
                index,
                row[chr_col],
                start,
                end,
                row[protein_col],
                uniprot_col,
                pos_column,
            )

        output_path = config["unique_signal_table"]["file_path"].replace(".csv", f"{suffix}.csv")
        unique_signal_table.to_csv(output_path, index=False, sep=sep)

    conn.close()
    print(f"Results updated and saved to {output_path}")


def process_table(conn, table, start_col, end_col, protein_col, uniprot_col, pos_column):
    for index, row in table.iterrows():
        process_table_row(
            conn,
            table,
            index,
            row["chr"],
            row[start_col],
            row[end_col],
            row[protein_col],
            row[uniprot_col],
            pos_column,
        )


def process_table_row(conn, table, index, chr, start, end, protein, uniprot, pos_column):
    # Step 1: Check for a SeqID match
    query_seqid = f"""
    SELECT {pos_column}, COHORT, pqtl_number_id FROM pQTL
    WHERE
        (chr = ? AND {pos_column} BETWEEN ? AND ? AND SeqID = ?)
    """
    cursor_seqid = conn.execute(query_seqid, (chr, start, end, protein))
    results_seqid = cursor_seqid.fetchall()

    if results_seqid:
        # If SeqID match is found
        pos_signals = [f"{chr}:{res[0]}" for res in results_seqid]
        cohorts = [str(res[1]) for res in results_seqid]
        number_ids = [str(res[2]) for res in results_seqid]

        table.at[index, "soma_match"] = "YES"
        table.at[index, "soma_matching_study"] = ",".join(cohorts)
        table.at[index, "soma_matching_number_ids"] = ",".join(number_ids)
        table.at[index, "soma_matching_signals"] = ",".join(pos_signals)

    query_uniprot = f"""
    SELECT {pos_column}, COHORT, pqtl_number_id FROM pQTL
    WHERE
        (chr = ? AND {pos_column} BETWEEN ? AND ? AND UniProt = ?)
    """
    cursor_uniprot = conn.execute(query_uniprot, (chr, start, end, uniprot))
    results_uniprot = cursor_uniprot.fetchall()

    if results_uniprot:
        # If UniProt match is found
        pos_signals = [f"{chr}:{res[0]}" for res in results_uniprot]
        cohorts = [str(res[1]) for res in results_uniprot]
        number_ids = [str(res[2]) for res in results_uniprot]

        table.at[index, "uniprot_match"] = "YES"
        table.at[index, "unip_matching_study"] = ",".join(cohorts)
        table.at[index, "unip_matching_number_ids"] = ",".join(number_ids)
        table.at[index, "unip_matching_signals"] = ",".join(pos_signals)


if __name__ == "__main__":
    backward_literature_review_annotation()
