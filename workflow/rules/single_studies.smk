rule appending_single_studies_results:
    input:
        rules.backward_literature_LB.output,
    output:
        ws_path("mapped_annotated_LB_lit_annotated_single_studies.csv"),
    conda:
        "../envs/single_studies.yml"
    resources:
        runtime=lambda wc, attempt: attempt * 60,
    shell:
        """workflow/scripts/single_studies/appending.sh"""
