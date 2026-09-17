# ----------------------------------------------------------------------
#  Project handling
# ----------------------------------------------------------------------
# Default project (used if no .project file exists)
DEFAULT_PROJECT ?= example

# .project file – contains the name of the active project
PROJECT_FILE := .project

# Resolve the active project:
#   1️⃣ If the file exists, read its contents.
#   2️⃣ Otherwise fall back to DEFAULT_PROJECT.
PROJECT := $(or $(shell cat $(PROJECT_FILE) 2>/dev/null),$(DEFAULT_PROJECT))

# Full config‑file path derived from the project name
CONFIGFILE := config/config.$(PROJECT).yaml

# ----------------------------------------------------------------------
#  Common Snakemake command parts
# ----------------------------------------------------------------------
CONDA_ENV_DIR=$(shell dirname ${CONDA_EXE})
HN=$(shell hostname | sed "s/[0-9]//g")
ifeq ($(HN),$(filter $(HN),cnode gnode hnode))
  CONDA_ENV_NAME=/exchange/healthds/software/envs/snakemake
else
  CONDA_ENV_NAME=snakemake
endif

# Targets list (unchanged)
TARGETS=dependencies dag run unlock

# ----------------------------------------------------------------------
#  Top level
# ----------------------------------------------------------------------
all:
	@echo "Try one of: ${TARGETS}"
	@echo "Current project: $(PROJECT) (config → $(CONFIGFILE))"

# ----------------------------------------------------------------------
#  DAG
# ----------------------------------------------------------------------
dag:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	snakemake --dag | dot -Tsvg > dag.svg

# ----------------------------------------------------------------------
#  Environment
# ----------------------------------------------------------------------
dependencies:
	mamba env update -n snakemake --file environment.yml

dev-dependencies: dependencies
	mamba env update -n snakemake --file environment_dev.yml

# ----------------------------------------------------------------------
#  Dry‑run
# ----------------------------------------------------------------------
dry-run:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	snakemake --sdm conda --dry-run --profile slurm --snakefile workflow/Snakefile

# ----------------------------------------------------------------------
#  Run
# ----------------------------------------------------------------------
pre-commit:
	if [ ! -f .git/hooks/pre-commit ]; then pre-commit install; fi
	pre-commit run --all-files

local-run:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	snakemake --printshellcmds --sdm conda --sdm apptainer --cores 4 --snakefile workflow/Snakefile

run:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	snakemake --profile slurm --snakefile workflow/Snakefile

rerun:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	snakemake --profile slurm --snakefile workflow/Snakefile --rerun-incomplete --rerun-trigger mtime

unlock:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	snakemake --unlock

dockerfile_:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	snakemake --containerize --snakefile workflow/Snakefile > Dockerfile

# ----------------------------------------------------------------------
#  Project helpers – create/override the .project file
# ----------------------------------------------------------------------
project-metaanalysis:
	@echo "metaanalysis" > $(PROJECT_FILE)
	@echo "✅ Project set to 'metaanalysis' "

project-believe:
	@echo "believe" > $(PROJECT_FILE)
	@echo "✅ Project set to 'believe' "

project-%:
	@echo "$*" > $(PROJECT_FILE)
	@echo "✅ Project set to '$*' "

# ----------------------------------------------------------------------
#  Clean up helper (optional)
# ----------------------------------------------------------------------
clean-project:
	@rm -f $(PROJECT_FILE)
	@echo "🗑️  Removed $(PROJECT_FILE); next run will fall back to DEFAULT_PROJECT='$(DEFAULT_PROJECT)'"