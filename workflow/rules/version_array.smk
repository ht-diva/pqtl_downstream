rule version_array:
    input:
        rules.select_best_SNP_from_LocusBreaker.output.mapped,
    output:
        annotated=ws_path("mapped_annotated_LB.csv"),
    resources:
        runtimes=lambda wc, attempt: attempt * 20,
    params:
        path_to_targets_list=config.get("array_list_path"),
    conda:
        "../envs/r_environment.yml"
    shell:
        """
         sumstats_path=$(dirname $(dirname $(head -n 1 {params.sumstats_list})));
         Rscript workflow/scripts/version_mapping/s01_annotate_version.R \
            --input {input} \
            --array_path {params.path_to_targets_list} \
            --annot_output {output.annotated}
   """
