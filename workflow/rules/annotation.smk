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
         Rscript workflow/scripts/version_mapping/s01_annotate_version.R \
            --input {input} \
            --array_path {params.path_to_targets_list} \
            --annot_output {output.annotated}
   """


rule gene_and_protein_annotation:
    input:
        rules.version_array.output.annotated,
    output:
        annotated=ws_path("mapped_annotated_LB_gp_ann.csv"),
    resources:
        runtimes=lambda wc, attempt: attempt * 20,
    params:
        collected_loci=config.get("LB_file"),
        mapping_file=config.get("mapping_filepath"),
        gtf_file=config.get("gtf_file"),
        uniprot_file=config.get("uniprot_file"),
    conda:
        "../envs/r_environment.yml"
    shell:
        """
         Rscript workflow/scripts/gp_annotation/gene_annotation_.R \
            --input {input} \
            --lb_file {params.collected_loci} \
            --mapping_file {params.mapping_file} \
            --gtf_file {params.gtf_file} \
            --uniprot_file {params.uniprot_file} \
            --output {output.annotated}
   """


rule backward_literature_LB:
    input:
        rules.gene_and_protein_annotation.output.annotated,
    output:
        ws_path("mapped_annotated_LB_gp_ann_lit_annotated.csv"),
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
        ws_path("mapped_annotated_LB_gp_ann_lit_annotated.hotspot_ann.csv"),
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
        "--save_results {params.save_results} "


rule appending_single_studies_results:
    input:
        rules.hostspot_finder.output,
    output:
        ws_path(
            "mapped_annotated_LB_gp_ann_lit_annotated.hotspot_ann_single_studies_{single_studies}.csv"
        ),
    conda:
        "../envs/single_studies.yml"
    log:
        ws_path("logs/single_studies_{single_studies}.log"),
    resources:
        runtime=lambda wc, attempt: attempt * 60,
    script:
        """../scripts/single_studies/appending.sh"""
