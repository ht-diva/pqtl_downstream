rule select_best_SNP_from_LocusBreaker:
    input:
        lb=rules.collect_loci.output.ofile,
        mapping=config.get("mapping_filepath"),
    output:
        MR=ws_path("MR_instruments_best_snps_from_LB.txt"),
        mapped=ws_path("mapped_LB.txt"),
        annotated=ws_path("mapped_annotated_LB.txt"),
    resources:
        runtimes=lambda wc, attempt: attempt * 20,
    params:
        NEF=config.get("params").get("nef"),
        sumstats_path=config.get("sum_stat_path"),
        path_to_targets_list=config.get("array_list_path"),
    conda:
        "../envs/r_environment.yml"
    shell:
        """
         Rscript workflow/scripts/s01_best_snp_locus_breaker_for_MR.R \
            --input {input.lb} \
            --path {params.sumstats_path} \
            --array_path {params.path_to_targets_list} \
            --mapping {input.mapping} \
            --NEF {params.NEF} \
            --map_output {output.mapped} \
            --MR_output {output.MR} \
            --annot_output {output.annotated}
   """
