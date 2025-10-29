# ==========================================================
# Paper A — Reproducible SDM Pipeline (H2O vs SSDM)
# Makefile orchestrating all steps
# ==========================================================

# Default parameters (can override: e.g. `make RUN=A Q=0.75`)
RUN ?= B
MODE ?= REPRO
Q ?= 0.75

# Script paths
SCRIPTS = scripts
RSCRIPT = Rscript

# ----------------------------------------------------------
# Environment setup (only once per system)
# ----------------------------------------------------------
setup:
	@echo ">>> Restoring R environment from renv.lock"
	@$(RSCRIPT) -e "install.packages('renv', repos='https://cloud.r-project.org'); renv::activate(); renv::restore()"

# ----------------------------------------------------------
# DBSCAN thinning + summary (optional upstream steps)
# ----------------------------------------------------------
thin:
	@$(RSCRIPT) $(SCRIPTS)/02_dbscan_thin_degrees.R --mode $(MODE)
	@$(RSCRIPT) $(SCRIPTS)/02a_summarize_thin.R

# ----------------------------------------------------------
# H2O and SSDM training (replicate workflows)
# ----------------------------------------------------------
rep-h2o:
	@echo ">>> Running H2O replicates for run $(RUN)"
	@$(RSCRIPT) $(SCRIPTS)/03a_h2o_replicates.R --run $(RUN) --mode $(MODE)

rep-ssdm:
	@echo ">>> Running SSDM replicates for run $(RUN)"
	@$(RSCRIPT) $(SCRIPTS)/04a_ssdm_replicates.R --run $(RUN) --mode $(MODE)

# ----------------------------------------------------------
# Between-method comparison from replicate results
# ----------------------------------------------------------
compare:
	@echo ">>> Comparing H2O vs SSDM for run $(RUN) (q=$(Q))"
	@$(RSCRIPT) $(SCRIPTS)/05a_h2o_vs_ssdm_results.R --run $(RUN) --mode $(MODE) --q $(Q)

# ----------------------------------------------------------
# Figures (no recomputation)
# ----------------------------------------------------------
paper:
	@echo ">>> Generating paper-ready figures from paper_results/"
	@$(RSCRIPT) $(SCRIPTS)/06_uncertainity.R

# ----------------------------------------------------------
# Full pipeline (replicates + comparison + paper)
# ----------------------------------------------------------
all: rep-h2o rep-ssdm compare paper

# ----------------------------------------------------------
# Clean up results
# ----------------------------------------------------------
clean:
	rm -rf results/* paper_results/* logs/*
	@echo ">>> Cleaned outputs."
