import pandas as pd
import sqlite3
import json
import click



def load_config(config_path):
    with open(config_path, 'r') as file:
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
def backward_literature_review_annotation(config_path, file_path, input_type):
    """
    Main function to orchestrate the literature review annotation process. This function
    handles the connection to the database, determines which type of data table (regional associations or unique signal)
    to process based on input, and ensures that the database connection is closed after processing.

    Inputs:
    - config_path: str
        Path to the JSON configuration file which contains all the necessary settings such as database path,
        reference genome version, and details about data tables (file paths, column names, etc.).
    - input_type: str
        Specifies the type of table to process. This can be either 'regional_table' or 'unique_signal_table'.
        The choice determines which dataset is loaded and processed, each requiring specific configurations.

    Outputs:
    - No direct return value, but this function triggers data processing which results in writing
      updated results to an updated CSV file. It updates tables
      with literature matching results and potentially secondary lookup results if no direct matches are found.

    """
    # Load the configuration file
    config = load_config(config_path)
    # overwrite the file_path
    config[input_type]['file_path'] = file_path

    # Connect to the SQLite database
    conn = sqlite3.connect(config["database_path"])

    # Determine the position column based on the reference genome version
    reference_genome = config["reference_genome"]
    pos_column = "pos37" if reference_genome == "37" else "pos38"

    # Process data based on the input type: regional table or unique signal table
    if input_type == 'regional_table':
        process_regional_table(conn, config, pos_column)
    elif input_type == 'unique_signal_table':
        process_unique_signal_table(conn, config, pos_column)

    # Close the database connection after processing
    conn.close()



def process_regional_table(conn, config, pos_column):
    regional_table = pd.read_csv(config["regional_table"]["file_path"])
    regional_table["matched"] = "NO"
    regional_table["matching_signals"] = ""
    regional_table["matching_study"] = ""
    regional_table["matching_number_ids"] = ""
    regional_table["secondary_found"] = "NO"
    regional_table["secondary_associated_targets"] = ""

    cohorts = config["regional_table"].get("cohorts", "all")
    start_col = config["regional_table"]["start_column"]
    end_col = config["regional_table"]["end_column"]
    protein_col = config["regional_table"]["protein_column"]

    process_table(conn, regional_table, start_col, end_col, protein_col, pos_column, cohorts)

    output_path = config["regional_table"]["file_path"].replace('.csv', '_lit_annotated.csv')
    regional_table.to_csv(output_path, index=False)
    print(f"Results updated and saved to {output_path}")


def process_unique_signal_table(conn, config, pos_column):
    unique_signal_table = pd.read_csv(config["unique_signal_table"]["file_path"])
    unique_signal_table["matched"] = "NO"
    unique_signal_table["matching_signals"] = ""
    unique_signal_table["matching_study"] = ""
    unique_signal_table["matching_number_ids"] = ""
    unique_signal_table["secondary_found"] = "NO"
    unique_signal_table["secondary_associated_targets"] = ""

    cohorts = config["unique_signal_table"].get("cohorts", "all")
    chr_col = config["unique_signal_table"]["chr_column"]
    pos_col = config["unique_signal_table"]["pos_column"]
    protein_col = config["unique_signal_table"]["protein_column"]
    region_size = config["unique_signal_table"]["region_size"]

    for index, row in unique_signal_table.iterrows():
        start = row[pos_col] - region_size
        end = row[pos_col] + region_size
        process_table_row(conn, unique_signal_table, index, row[chr_col], start, end, row[protein_col], pos_column, cohorts, region_size)

    output_path = config["unique_signal_table"]["file_path"].replace('.csv', '_lit_annotated.csv')
    unique_signal_table.to_csv(output_path, index=False)
    print(f"Results updated and saved to {output_path}")



def process_table(conn, table, start_col, end_col, protein_col, pos_column, cohorts):
    for index, row in table.iterrows():
        process_table_row(conn, table, index, row['chr'], row[start_col], row[end_col], row[protein_col], pos_column, cohorts, row[start_col] - row[end_col])

def process_table_row(conn, table, index, chr, start, end, protein, pos_column, cohorts, region_size):
    primary_found = execute_primary_lookup(conn, table, index, chr, start, end, protein, pos_column, cohorts)
    if not primary_found:
        execute_secondary_lookup(conn, table, index, chr, start, end, pos_column, region_size)


def execute_primary_lookup(conn, table, index, chr, start, end, protein, pos_column, cohorts):
    cohort_filter = f"AND COHORT IN ({', '.join('?' for _ in cohorts)})" if cohorts != "all" else ""
    query = f'''
    SELECT {pos_column}, COHORT, pqtl_number_id FROM pQTL
    WHERE
        (chr = ? AND {pos_column} BETWEEN ? AND ? AND (SeqID = ? OR OlinkID = ? OR UniProt = ?))
        {cohort_filter}
    '''
    query_params = [chr, start, end, protein, protein, protein] + (cohorts if cohorts != "all" else [])
    cursor = conn.execute(query, query_params)
    results = cursor.fetchall()

    if results:
        update_table_with_results(table, index, results, chr, pos_column)
        return True
    return False

def execute_secondary_lookup(conn, table, index, chr, start, end, pos_column, region_size):
    query = f'''
    SELECT {pos_column}, COHORT, SeqID FROM pQTL
    WHERE
        (chr = ? AND {pos_column} BETWEEN ? AND ?)
    '''
    cursor = conn.execute(query, (chr, start, end))
    results = cursor.fetchall()

    if results:
        targets = set(res[2] if res[2] is not None else 'Unknown' for res in results)
        table.at[index, "secondary_found"] = "YES"
        table.at[index, "secondary_associated_targets"] = ", ".join(map(str,targets))
    else:
        table.at[index, "secondary_found"] = "NO"
        table.at[index, "secondary_associated_targets"] = ""


def update_table_with_results(table, index, results, chr, pos_column):
    pos_signals = [f"{chr}:{res[0]}" for res in results]
    cohorts = [str(res[1]) for res in results]
    number_ids = [str(res[2]) for res in results]

    table.at[index, "matched"] = "YES"
    table.at[index, "matching_signals"] = ",".join(pos_signals)
    table.at[index, "matching_study"] = ",".join(cohorts)
    table.at[index, "matching_number_ids"] = ",".join(number_ids)

if __name__ == "__main__":
    backward_literature_review_annotation()
