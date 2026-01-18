# ==========================================================
# Paper A — Reproducible SDM Pipeline (H2O vs SSDM)
# Makefile orchestrating all steps
# ==========================================================

# Default parameters (can override: e.g. `make RUN=A Q=0.75`)
RUN ?= B
MODE ?= REPRO
Q ?= 0.75

# Script paths
SCRIPTS_HS = scripts_Habitat_Suitability
SCRIPTS_ME = scripts_Movement_Ecology
RSCRIPT = Rscript

# ----------------------------------------------------------
# Environment setup
# ----------------------------------------------------------
setup:
	@echo ">>> Restoring R environment from renv.lock"
	@$(RSCRIPT) -e "install.packages('renv', repos='https://cloud.r-project.org'); renv::activate(); renv::restore()"

# ----------------------------------------------------------
# DBSCAN thinning
# ----------------------------------------------------------
thin:
	@$(RSCRIPT) $(SCRIPTS_ME)/02_spatial_thinning.R --run $(RUN) --mode $(MODE)

# ----------------------------------------------------------
# H2O and SSDM training
# ----------------------------------------------------------
train-h2o:
	@echo ">>> Running H2O training for run $(RUN)"
	@$(RSCRIPT) $(SCRIPTS_HS)/03_h2o_train.R --run $(RUN) --mode $(MODE)

train-ssdm:
	@echo ">>> Running SSDM training for run $(RUN)"
	@$(RSCRIPT) $(SCRIPTS_HS)/04_ssdm_train.R --run $(RUN) --mode $(MODE)

# ----------------------------------------------------------
# Comparison & Appendix
# ----------------------------------------------------------
compare:
	@echo ">>> Comparing H2O vs SSDM for run $(RUN)"
	@$(RSCRIPT) $(SCRIPTS_HS)/05_h2o_vs_ssdm.R --run $(RUN) --mode $(MODE)

appendix:
	@echo ">>> Generating appendix figures"
	@$(RSCRIPT) $(SCRIPTS_HS)/07_appendix.R

# ----------------------------------------------------------
# Full pipeline
# ----------------------------------------------------------
all: thin train-h2o train-ssdm compare appendix

# ----------------------------------------------------------
# Clean up results
# ----------------------------------------------------------
clean:
	rm -rf results/H2O/* results/SSDM/* results/compare/* results/appendix/* logs/*
	@echo ">>> Cleaned outputs."
