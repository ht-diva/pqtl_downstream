FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="58307e12d822423826be0455d8d55caf05de14123b3e78fd1fd2782a8a759763"

# Step 1: Retrieve conda environments

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
#   prefix: /conda-envs/6b98aaef539fd8905de89427d31bbd7e
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
#     - r-R.utils=2.12.3
#     - bioconductor-iranges=2.36.0
RUN mkdir -p /conda-envs/6b98aaef539fd8905de89427d31bbd7e
COPY workflow/envs/r_environment.yml /conda-envs/6b98aaef539fd8905de89427d31bbd7e/environment.yaml

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

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/1bcfde7b7f5d72aa0ee14f0145284fff --file /conda-envs/1bcfde7b7f5d72aa0ee14f0145284fff/environment.yaml && \
    mamba env create --prefix /conda-envs/6b98aaef539fd8905de89427d31bbd7e --file /conda-envs/6b98aaef539fd8905de89427d31bbd7e/environment.yaml && \
    mamba env create --prefix /conda-envs/e68b6dcacdcf9733ffce40eff8569902 --file /conda-envs/e68b6dcacdcf9733ffce40eff8569902/environment.yaml && \
    mamba clean --all -y
