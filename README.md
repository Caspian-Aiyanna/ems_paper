# 🐘 BTEH – Elephant Habitat Modelling in restored landscapes

This repository contains the full, reproducible modelling pipeline of the BTEH project. It quantifies **habitat suitability, temporal change and comparing two SDM techniques** for African elephants before and after fence removal in Kariega Game Reserve, South Africa.

All workflows are implemented in **R** (with optional Google Earth Engine extraction scripts) and are designed to run both locally and on **HPC clusters** via SLURM or PBS job schedulers.

---

## 🌿 Project Overview

The pipeline integrates multiple modelling frameworks:

| Framework | Description |
|------------|--------------|
| **H2O AutoML** | Machine-learning ensemble (GBM, XGBoost, DNN, DRF, etc - upto 100 models) for habitat suitability |
| **SSDM** | Classical stacked species distribution modelling (GLM, GAM, RF, ANN, GBM) |
| **Comparison & Metrics** | Between-method (H2O vs SSDM), temporal change, and hotspot analyses |

All outputs are automatically written under `results/` and logged under `logs/`.

---

## 📁 Folder Structure

```
.
├── BTEH/
│   ├── config.yml                # Global config (mode, cores, etc.)
│   ├── Makefile                  # Optional automation targets
│   ├── scripts_Habitat_Suitability/ # Core SDM Pipeline (H2O, SSDM, Appendix)
│   ├── scripts_Movement_Ecology/    # Behavioral & Temporal analysis (Thinning, Panels, GSL)
│   ├── R/                           # Utility functions
│   ├── config.yml                   # Global config (mode, cores, etc.)
│   ├── Makefile                     # Automation targets
│   ├── data/
│   │   ├── clean/       # Pre-processed telemetry datasets
│   │   ├── envi/        # Environmental sets
│   │   ├── occ/         # Occurrence replicates
│   │   └── shp/         # Boundary shapefiles
│   ├── results/          # Generated map outputs
│   ├── logs/             # Execution logs
│   ├── plans/            # Variable selection plans
├── renv.lock                 # Environment lockfile
└── results/                  # Final manuscript assets
```

---

## ⚙️ Environment Setup

1. **Clone the repository**
   ```bash
   git clone https://github.com/Caspian-Aiyanna/ems_paper.git
   cd ems_paper
   ```

2. **Restore the R environment**
   ```bash
   Rscript -e "renv::restore(prompt = FALSE)"
   ```

3. **(Optional)** Check directory structure
   ```r
   fs::dir_tree('.', recurse = 2)
   ```

---

## 🚀 Running the Pipeline

### 🧩 Option 1 — Reproducible (“REPRO”) Run

Use **single-core deterministic** execution to eliminate randomness.  
Ensures rigorous and reproducible results.

```bash
sbatch hpc/BTEH_REPRO.slurm
```

Or locally:

```bash
Rscript scripts_Habitat_Suitability/03_h2o_train.R --run A --mode REPRO
Rscript scripts_Habitat_Suitability/04_ssdm_train.R --run A --mode REPRO
Rscript scripts_Habitat_Suitability/05_h2o_vs_ssdm.R --run A --mode REPRO

Rscript scripts_Habitat_Suitability/03_h2o_train.R --run B --mode REPRO
Rscript scripts_Habitat_Suitability/04_ssdm_train.R --run B --mode REPRO
Rscript scripts_Habitat_Suitability/05_h2o_vs_ssdm.R --run B --mode REPRO
```

**Outputs:**  
- `results/H2O/<RUN>/<SP>/prediction_<SP>.tif`  
- `results/SSDM/<RUN>/<SP>/ESDM_<SP>.tif`  
- `results/compare/h2o_vs_ssdm/<RUN>/metrics/…`  
- Logs under `logs/03_*.log`, `logs/04_*.log`, `logs/05_*.log`

---

### ⚡ Option 2 — Fast (“FAST”) Run

Parallelized execution for development and testing.  
Results are near-identical but may vary slightly due to parallel randomness.

```bash
sbatch hpc/BTEH_FAST.slurm
```

Or locally (multi-core machine):

```bash
Rscript scripts_Habitat_Suitability/03_h2o_train.R --run A --mode FAST
Rscript scripts_Habitat_Suitability/04_ssdm_train.R --run A --mode FAST
Rscript scripts_Habitat_Suitability/05_h2o_vs_ssdm.R --run A --mode FAST
```

---

### 🧬 Option 3 — Species-Array Mode (HPC only)

Runs each elephant (E1B–E6B) as a separate array job:

```bash
sbatch hpc/BTEH_FAST_array.slurm
```

This distributes species across nodes and merges results automatically.

---

## 🧾 Log Files

All scripts log progress and warnings to the `logs/` folder:
- `03_h2o_train_A.log`, `04_ssdm_train_B.log`, etc.
- `05_compare_methods_A.log` and `05_compare_methods_B.log`
- Environment information and timing details are captured for reproducibility.

---

## 📊 Expected Outputs

- `results/H2O/` – AutoML rasters, models, leaderboards  
- `results/SSDM/` – Ensemble rasters, algorithm summaries  
- `results/compare/` – Metrics, hotspot overlaps, and temporal change maps  
- `results/appendix/` – Supplementary figures (base learner weights, etc.)  
- `results/panels/` – Multi-pane composite panels  

Each result folder includes `.csv` summaries and `.tif` rasters ready for visualization.

---

## 🧠 Notes for HPC Users

- **Modules:** Load `R/4.3.2`, `GDAL`, `GEOS`, and `PROJ` if required.  
- **Logs:** Check progress with `tail -f logs/BTEH_REPRO_*.out`.  
- **Monitoring:** Use `squeue -u $USER` (SLURM) or `qstat` (PBS).  
- **Storage:** Prefer `/scratch` for heavy intermediate files.

---

## 📚 Citation

If you use this workflow or data structure, please cite:

> Harin Aiyanna C R; Francesco Pirotti; Brooke Friswold; Antoinette van de Water. (2025). *Movement Meets Machine Learning: Dual framework testing for Predicting Elephant Habitat Suitability*.  
> Movement Ecology. (In preparation)

---

## 🧰 Contact

**Author:** Harin Aiyanna C R
**Institution:** [Interdepartmental Research Center of Geomatics (CIRGEO), University of Padova and Bring The Elephant Home]  
**Email:** [harinaiyanna.cheriyandaraveendra@phd.unipd.it]  
**GitHub:** [https://github.com/Caspian-Aiyanna](https://github.com/Caspian-Aiyanna)

---

### 🏁 Quick Summary

| Mode | Purpose | HPC script | CPU | Typical runtime |
|------|----------|-------------|-----|-----------------|
| **REPRO** | Final, deterministic results | `hpc/BTEH_REPRO.slurm` | 1 | 24–48 h |
| **FAST** | Development, parallel | `hpc/BTEH_FAST.slurm` | 8–16 | 3–8 h |
| **FAST-ARRAY** | Multi-species parallel | `hpc/BTEH_FAST_array.slurm` | 8 each | 1–4 h |

---

*This repository embodies the principles of open, reproducible ecological modelling and provides a modular foundation for cross-framework SDM benchmarking.*
