rule collapsing:
    input:
        rules.backward_literature_LB.output,
    output:
        collapsed=ws_path(mapped_LB_collapsed),
    resources:
        runtime=lambda wc, attempt: attempt * 20,
    params:
        mapping_file=config.get("mapping_filepath"),
    conda:
        "../envs/r_environment.yml"
    shell:
        """
         Rscript workflow/scripts/collapsing/s01_collapsing.R \
            --input {input} \
            --mapping {params.mapping_file} \
            --output {output.collapsed}
   """
