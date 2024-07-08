from pathlib import Path
import pandas as pd


# Define input for the rules
data = []
with open(config["sumstats_list"], "r") as fp:
    lines = fp.readlines()

for line in lines:
    p = Path(line.strip())
    seqid = ".".join(p.stem.split(".")[:3])  # seqid = p.stem.split("_")[3]
    data.append((seqid, str(p)))

analytes = (
    pd.DataFrame.from_records(data, columns=["seqid", "sumstat_path"])
    .set_index("seqid", drop=False)
    .sort_index()
)


def get_final_output():
    final_output = []

    if config.get("input") == "LB":
        final_output.append(ws_path("break/collected_loci_excluding_mhc.csv")),
        final_output.append(ws_path("MR_instruments_best_snps_from_LB.txt")),
        final_output.append(ws_path("mapped_LB.txt")),
        final_output.append(ws_path("mapped_annotated_LB.txt")),

    return final_output


def get_sumstats(wildcards):
    return analytes.loc[wildcards.seqid, "sumstat_path"]


# define the functions generating files' path
def ws_path(file_path):
    return str(Path(config.get("workspace_path"), file_path))
