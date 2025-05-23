if config.get("input") == "LB":
    input_MR = config.get("LB_file")
elif config.get("input") == "run_LB":
    input_MR = rules.collect_loci.output.ofile


rule select_best_SNP_from_LocusBreaker:
    input:
        lb=input_MR,
        mapping=config.get("mapping_filepath"),
    output:
        MR=ws_path("MR_instruments_best_snps_from_LB.txt"),
        mapped=ws_path(select_best_SNP_from_LocusBreaker_output),
    resources:
        runtime=lambda wc, attempt: attempt * 20,
    params:
        NEF=config.get("params").get("nef"),
        sumstats_list=config.get("sumstats_list"),
    conda:
        "../envs/r_environment.yml"
    shell:
        """
         sumstats_path=$(dirname $(dirname $(head -n 1 {params.sumstats_list})));
         Rscript workflow/scripts/MR/s01_best_snp_locus_breaker_for_MR.R \
            --input {input.lb} \
            --path $sumstats_path \
            --mapping {input.mapping} \
            --NEF {params.NEF} \
            --map_output {output.mapped} \
            --MR_output {output.MR}
   """
