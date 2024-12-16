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


select_best_SNP_from_LocusBreaker_output = "mapped_LB.csv"

outputs = {
    "gp": "_gp_ann",
    "va": "_va_ann",
    "bl": "_bl_ann",
    "hf": "_hf_ann",
    "as": "_as_ann",
}

annotation_output = select_best_SNP_from_LocusBreaker_output
annotation_outputs = {}
for suffix, output in outputs.items():
    annotation_output = annotation_output.replace(".csv", f"{output}.csv")
    annotation_outputs[suffix] = annotation_output


def get_final_output():
    final_output = []

    if config.get("input") == "run_LB":
        final_output.append(rules.collect_loci.output.ofile),
        final_output.append(rules.select_best_SNP_from_LocusBreaker.output.MR),
        final_output.append(rules.select_best_SNP_from_LocusBreaker.output.mapped),
        final_output.append(rules.heterogenous_filter.output),
        final_output.append(rules.collapse.output),
        final_output.extend(
            expand(
                rules.appending_single_studies_results.output,
                single_studies=config.get("single_studies"),
            )
        ),

    if config.get("input") == "LB":
        final_output.append(rules.select_best_SNP_from_LocusBreaker.output.MR),
        final_output.append(rules.select_best_SNP_from_LocusBreaker.output.mapped),
        final_output.append(rules.heterogenous_filter.output),
        final_output.append(rules.collapse.output),
        final_output.extend(
            expand(
                rules.appending_single_studies_results.output,
                single_studies=config.get("single_studies"),
            )
        ),

    return final_output


def get_sumstats(wildcards):
    return analytes.loc[wildcards.seqid, "sumstat_path"]


# define the functions generating files' path
def ws_path(file_path):
    return str(Path(config.get("workspace_path"), file_path))
