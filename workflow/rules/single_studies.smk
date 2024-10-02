rule appending_single_studies_results:
    input:
        rules.backward_literature_LB.output,
    output:
        expand(
            ws_path(
                "mapped_annotated_LB_lit_annotated_single_studies_{single_studies}.csv",
                single_studies=config.get("single_studies"),
            )
        ),
    conda:
        "../envs/single_studies.yml"
    log:
        ws_path("logs/single_studies.log"),
    resources:
        runtime=lambda wc, attempt: attempt * 60,
    script:
        """../scripts/single_studies/appending.sh"""
