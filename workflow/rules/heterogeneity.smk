rule heterogenous_filter:
    input:
        rules.select_best_SNP_from_LocusBreaker.output.annotated,
    output:
        ws_path("heterogenous_LB.txt"),
    resources:
        runtimes=lambda wc, attempt: attempt * 20,
    params:
        NEF=config.get("params").get("nef"),
        Isquare=config.get("params").get("Isquare"),
    conda:
        "../envs/r_environment.yml"
    shell:
        """
         Rscript workflow/scripts/heterogeneity/s03_heterogenous_filter.R \
            --input {input} \
            --NEF {params.NEF} \
            --het_output {output} \
            --Isquare_thresh {params.Isquare}
   """
