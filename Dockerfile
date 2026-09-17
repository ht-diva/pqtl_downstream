FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="f431b0eee4992d93e3e14ebefe6b1187b03cec80fb904ba2eabe98da1503d695"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/locus_breaker.yml
#   prefix: /conda-envs/1bcfde7b7f5d72aa0ee14f0145284fff
#   name: locus_breaker
#   channels:
#     - conda-forge
#     - defaults
#     - R
#   dependencies: #R>=4.3
#     - r-base=4.3.3
#     - r-optparse=1.7.4
#     - r-tidyverse=2.0.0
#     - r-data.table=1.15.2
#     - r-R.utils=2.12.3
RUN mkdir -p /conda-envs/1bcfde7b7f5d72aa0ee14f0145284fff
COPY workflow/envs/locus_breaker.yml /conda-envs/1bcfde7b7f5d72aa0ee14f0145284fff/environment.yaml

# Conda environment:
#   source: workflow/envs/r_environment.yml
#   prefix: /conda-envs/43aa15f1edfb65818efdbb029ca8ad66
#   name: r_MR
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#     - R
#   dependencies: #R>=4.3
#     - r-base=4.3.3
#     - r-optparse=1.7.4
#     - r-tidyverse=2.0.0
#     - r-dplyr=1.1.4
#     - r-data.table=1.15.4
#     - r-R.utils=2.12.3
#     - bioconductor-iranges=2.36.0
#     - bioconductor-genomicranges=1.54.1
#     - bioconductor-rtracklayer=1.62.0
RUN mkdir -p /conda-envs/43aa15f1edfb65818efdbb029ca8ad66
COPY workflow/envs/r_environment.yml /conda-envs/43aa15f1edfb65818efdbb029ca8ad66/environment.yaml

# Conda environment:
#   source: workflow/envs/single_studies.yml
#   prefix: /conda-envs/14112df28a285a4ce6f6e5a42f114900
#   name: htslib
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - htslib=1.20
RUN mkdir -p /conda-envs/14112df28a285a4ce6f6e5a42f114900
COPY workflow/envs/single_studies.yml /conda-envs/14112df28a285a4ce6f6e5a42f114900/environment.yaml

# Conda environment:
#   source: workflow/scripts/backward_literature/environment.yml
#   prefix: /conda-envs/e68b6dcacdcf9733ffce40eff8569902
#   name: backward_literature
#   channels:
#     - conda-forge
#     - defaults
#   dependencies:
#     - python<=3.10
#     - pip==24
#     - pip:
#         - pandas==2.2
#         - pyarrow==15.0
#         - click==8.1.7
RUN mkdir -p /conda-envs/e68b6dcacdcf9733ffce40eff8569902
COPY workflow/scripts/backward_literature/environment.yml /conda-envs/e68b6dcacdcf9733ffce40eff8569902/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/1bcfde7b7f5d72aa0ee14f0145284fff --file /conda-envs/1bcfde7b7f5d72aa0ee14f0145284fff/environment.yaml && \
    conda env create --prefix /conda-envs/43aa15f1edfb65818efdbb029ca8ad66 --file /conda-envs/43aa15f1edfb65818efdbb029ca8ad66/environment.yaml && \
    conda env create --prefix /conda-envs/14112df28a285a4ce6f6e5a42f114900 --file /conda-envs/14112df28a285a4ce6f6e5a42f114900/environment.yaml && \
    conda env create --prefix /conda-envs/e68b6dcacdcf9733ffce40eff8569902 --file /conda-envs/e68b6dcacdcf9733ffce40eff8569902/environment.yaml && \
    conda clean --all -y
