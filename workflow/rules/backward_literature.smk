rule backward_literature_LB:
    input:
        rules.version_array.output.annotated,
    output:
        ws_path("mapped_annotated_LB_lit_annotated.csv"),
    conda:
        "../scripts/backward_literature/environment.yml"
    params:
        input_type="regional_table",
        config_file=config.get("BL_config_file"),
    resources:
        runtime=lambda wc, attempt: attempt * 60,
    shell:
        "python workflow/scripts/backward_literature/Backward_Literature.py "
        "-c {params.config_file} "
        "-f {input} "
        "-t {params.input_type} "
