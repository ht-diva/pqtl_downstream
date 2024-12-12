rule gene_and_protein_annotation:
    input:
        rules.select_best_SNP_from_LocusBreaker.output.mapped,
    output:
        annotated=ws_path(annotation_outputs["gp"]),
    resources:
        runtime=lambda wc, attempt: attempt * 60,
    params:
        mapping_file=config.get("mapping_filepath"),
        gtf_file=config.get("gtf_file"),
    conda:
        "../envs/r_environment.yml"
    shell:
        """
         Rscript workflow/scripts/gp_annotation/gene_annotations.R \
            --input {input} \
            --mapping {params.mapping_file} \
            --gtf_file {params.gtf_file} \
            --output {output.annotated}
   """


rule version_array:
    input:
        rules.gene_and_protein_annotation.output.annotated,
    output:
        annotated=ws_path(annotation_outputs["va"]),
    resources:
        runtime=lambda wc, attempt: attempt * 20,
    params:
        path_to_targets_list=config.get("array_list_path"),
    conda:
        "../envs/r_environment.yml"
    shell:
        """
         Rscript workflow/scripts/version_mapping/s01_annotate_version.R \
            --input {input} \
            --array_path {params.path_to_targets_list} \
            --annot_output {output.annotated}
   """


rule backward_literature_LB:
    input:
        rules.version_array.output.annotated,
    output:
        ws_path(annotation_outputs["bl"]),
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


rule hostspot_finder:
    input:
        rules.backward_literature_LB.output,
    output:
        ws_path(annotation_outputs["hf"]),
    conda:
        "../scripts/backward_literature/environment.yml"
    params:
        hotspot_window_size=config.get("hotspot_finder").get("hotspot_window_size"),
        chr_col=config.get("hotspot_finder").get("chr_col"),
        start_col=config.get("hotspot_finder").get("start_col"),
        end_col=config.get("hotspot_finder").get("end_col"),
        hotspot_threshold=config.get("hotspot_finder").get("hotspot_threshold"),
        lonespot_window_size=config.get("hotspot_finder").get("lonespot_window_size"),
        lonespot_threshold=config.get("hotspot_finder").get("lonespot_threshold"),
        save_results=config.get("hotspot_finder").get("save_results"),
    resources:
        runtime=lambda wc, attempt: attempt * 60,
    shell:
        "python workflow/scripts/hotspot/hotspot_finder.py "
        "--file_path {input} "
        "--hotspot_window_size {params.hotspot_window_size} "
        "--chr_col {params.chr_col} "
        "--start_col {params.start_col} "
        "--end_col {params.end_col} "
        "--hotspot_threshold {params.hotspot_threshold} "
        "--lonespot_window_size {params.lonespot_window_size} "
        "--lonespot_threshold {params.lonespot_threshold} "
        "--save_results "


rule appending_single_studies_results:
    input:
        rules.hostspot_finder.output,
    output:
        ws_path(annotation_outputs["as"].replace(".csv", "_{single_studies}.csv")),
    conda:
        "../envs/single_studies.yml"
    params:
        study="{single_studies}",
    log:
        ws_path("logs/single_studies_{single_studies}.log"),
    resources:
        runtime=lambda wc, attempt: attempt * 60,
    script:
        """../scripts/single_studies/appending.sh"""
