import click
import math
import pandas as pd


def hotspot_finder(file_path, hotspot_window_size, chr_col, start_col, end_col,
                   hotspot_threshold=50, lonespot_window_size=None, lonespot_threshold=None, save_results=False, separator=';'):
    hotspot_window_size = int(hotspot_window_size)
    hotspot_threshold = int(hotspot_threshold)
    lonespot_window_size = int(lonespot_window_size)
    lonespot_threshold = int(lonespot_threshold)

    df = pd.read_csv(file_path, sep=separator)

    df['hotspot_windows'] = df.apply(lambda row: calculate_covered_windows(row, hotspot_window_size, start_col, end_col), axis=1)

    # Now we create a dictionary to store the count of how many times each window is covered for each chromosome (hotspots)
    hotspot_count = {}

    # Iterate through each chromosome
    for chr_num in df[chr_col].unique():
        chr_df = df[df[chr_col] == chr_num]

        region_count = {}

        for covered_windows in chr_df['hotspot_windows']:
            for window in covered_windows:
                if window not in region_count:
                    region_count[window] = 0
                region_count[window] += 1

        # Store the region count for this chromosome
        hotspot_count[chr_num] = region_count

    # Identify the hotspots based on the hotspot_threshold
    df['hotspot'] = df.apply(lambda row: is_in_hotspot(row, hotspot_window_size, start_col, end_col, hotspot_count, hotspot_threshold), axis=1)

    # If lonespot_window_size and lonespot_threshold are provided, calculate lonespots
    if lonespot_window_size and lonespot_threshold:

        df['lonespot_windows'] = df.apply(lambda row: calculate_covered_windows(row, lonespot_window_size, start_col, end_col), axis=1)

        lonespot_count = {}

        for chr_num in df[chr_col].unique():

            chr_df = df[df[chr_col] == chr_num]
            chr_df = chr_df[chr_df['cis_or_trans'] == 'trans']
            region_count = {}


            for lonespot_windows in chr_df['lonespot_windows']:
                for window in lonespot_windows:
                    if window not in region_count:
                        region_count[window] = 0
                    region_count[window] += 1

            # Store the region count for this chromosome
            lonespot_count[chr_num] = region_count

        # Identify the lonespots based on the lonespot_threshold
        df['lonespot'] = df.apply(lambda row: is_in_lonespot(row, lonespot_window_size, start_col, end_col, lonespot_count, lonespot_threshold), axis=1)
        df['signal_count_lonespot'] = df.apply(lambda row: count_lonespot(row, lonespot_window_size, start_col, end_col, lonespot_count, lonespot_threshold), axis=1)
    if save_results:
        output_table_path = file_path.replace(".csv", "_hotspotted.csv")
        output_hotspot_dict_path = file_path.replace(".csv", "_hotspot_dict.csv")
        output_lonespot_dict_path = file_path.replace(".csv", "_lonespot_dict.csv")

        # Save the annotated table with hotspots and lonespots
        df.to_csv(output_table_path, sep=separator, index=False)
        print(f"Table with hotspots and lonespots saved to {output_table_path}")

        # Save the hotspot count dictionary
        dict_df = pd.DataFrame([(k, inner_k, v) for k, v_dict in hotspot_count.items() for inner_k, v in v_dict.items()],
                               columns=[chr_col, 'window', 'count'])
        dict_df.to_csv(output_hotspot_dict_path, sep="\t", index=False)
        print(f"Hotspot count dictionary saved to {output_hotspot_dict_path}")

        # Save the lonespot count dictionary
        if lonespot_window_size and lonespot_threshold:
            lonespot_dict_df = pd.DataFrame([(k, inner_k, v) for k, v_dict in lonespot_count.items() for inner_k, v in v_dict.items()],
                                            columns=[chr_col, 'window', 'count'])
            lonespot_dict_df.to_csv(output_lonespot_dict_path, sep="\t", index=False)
            print(f"Lonespot count dictionary saved to {output_lonespot_dict_path}")

    return df

def calculate_covered_windows(row, window_size, start_col, end_col):
    # Calculate the start_region and end_region based on window size
    start_region = math.ceil(row[start_col] / window_size)
    end_region = math.ceil(row[end_col] / window_size)

    # If start and end are in the same region, only one region is covered
    if start_region == end_region:
        return [start_region]

    # Otherwise, create a list of all regions between start and end (inclusive)
    return list(range(start_region, end_region + 1))

def is_in_hotspot(row, window_size, start_col, end_col, hotspot_count, hotspot_threshold):

    covered_windows = row['hotspot_windows']

    # Check if any of the covered windows in this row are a hotspot
    chr_num = row['chr']
    for window in covered_windows:
        if chr_num in hotspot_count and window in hotspot_count[chr_num]:
            if hotspot_count[chr_num][window] >= hotspot_threshold:
                return True
    return False

def is_in_lonespot(row, window_size, start_col, end_col, lonespot_count, lonespot_threshold):

    lonespot_windows = row['lonespot_windows']

    chr_num = row['chr']
    if row['cis_or_trans'] == 'trans':
        for window in lonespot_windows:
            if chr_num in lonespot_count and window in lonespot_count[chr_num]:
                if lonespot_count[chr_num][window] <= lonespot_threshold:
                    return True
    return False

def count_lonespot(row, window_size, start_col, end_col, lonespot_count, lonespot_threshold):

    lonespot_windows = row['lonespot_windows']

    chr_num = row['chr']
    count = 0
    for window in lonespot_windows:
        if chr_num in lonespot_count and window in lonespot_count[chr_num]:
            if lonespot_count[chr_num][window] > count:
                count = lonespot_count[chr_num][window]
    return count

@click.command()
@click.option("--file_path", required=True, help="Table file path")
@click.option("--hotspot_window_size", required=True, help="hotspot window size")
@click.option("--chr_col", required=True, default="chr", help="chromosome column name")
@click.option("--start_col", required=True, help="start column")
@click.option("--end_col", required=True, help="end column")
@click.option("--hotspot_threshold", default=50, required=True, help="hotspot threshold")
@click.option("--lonespot_window_size", default=None, required=True, help="lonespot window size")
@click.option("--lonespot_threshold", default=None, required=True, help="lonespot threshold")
@click.option("--save_results", required=True, default=False, is_flag=True, help="Table file path")
def main(file_path, hotspot_window_size, chr_col, start_col, end_col,
                   hotspot_threshold, lonespot_window_size, lonespot_threshold, save_results):

    hotspot_finder(
        file_path,
        hotspot_window_size, chr_col, start_col,
        end_col, hotspot_threshold, lonespot_window_size,
        lonespot_threshold, save_results)


if __name__ == "__main__":
    main()
