# Repository Overview

This repository contains the full pipeline for:
- building real‑data series (patents + papers),
- fitting parameters for the UMT‑TAP model,
- running simulations,
- comparing real data vs simulations,
- visualizing branching regimes.

## Directories
You will work with three directories, configured inside the notebooks:
- `data_dir_patents`: raw input data (you provide this).
- `output_dir`: derived outputs from the notebooks (auto‑generated).
- `model_dir`: model simulation outputs (auto‑generated). This is **`data_simulations/`**.

## Module 1 — Real Data Preparation

### Patents (`module_1_patents.ipynb`)
**Input (in `data_dir_patents`)**
- One CSV per year, named `YYYY.csv` (e.g., `1980.csv`, `1981.csv`, ...).
- Required columns:
  - `publication.date` (format `YYYYMMDD`)
  - `inventor` (list-like string of dicts, each with `name`)
  - `ipc` (list-like string)

**Processing**
- Excludes first authors with more than `max_patents_per_year` in any year.
- Builds cumulative series and frequency tables.

**Outputs (in `output_dir`)**
- `data_augmented_patents/` (per‑year CSVs with cumulative columns)
- `ipc_combos_freq.csv`
- `unique_authors_freq.csv`
- `unique_first_authors_freq.csv`
- `delta_rows_patents.csv`
- `reduced_patents_log.csv`
- `reduced_patents_lin.csv`
- `epo_author_dates_1980_2020.pkl`
- `epo_intervals.pkl`

### Papers (`module_1_papers.ipynb`)
**Input (in `data_dir_papers`)**
- OpenAlex subfield parquet files:
  - `merged_df_subfield_{ID}.parquet`
- Required columns:
  - `date`
  - `topics`
  - `keywords`
  - `authors` (first author used as proxy)

**Processing**
- Excludes first authors with more than `max_papers_per_year` in any year.
- Builds cumulative series and frequency tables.

**Outputs (in `output_dir`)**
- `delta_rows_papers.csv`
- `openalex_combo_freq.csv`
- `openalex_author_freq.csv`
- `reduced_openalex_log.csv`
- `reduced_openalex_lin.csv`
- `openalex_author_dates_1980_2020.pkl`
- `openalex_intervals.pkl` (if generated)

## Module 2 — Real Data Analysis & Fits
**Notebook:** `module_2_real_data.ipynb`

**Inputs**
- Reduced real‑data series from `output_dir`
- Author‑date pickles from Module 1
- Inter‑event interval files (if present)

**What it does**
- Builds the 2×2 fit panel (`D_w(D)` and `dt/dτ` vs `w`) for patents and papers.
- Computes yearly productivity metrics and plots rate figures.
- Computes and prints inter‑event time statistics (patents + papers).
- Builds author‑year activity heatmaps (1980–2020).

**Outputs**
- Figures only (no new data files required for later modules).

## Module 3 — UMT‑TAP Simulations (Micro–Macro)

### Compile
```bash
gcc -o model_micro_macro model_micro_macro.c -lm
```

### Run
```bash
./model_micro_macro T nu rho N0 a b w0
```

**Important:** use `T = 41` for both patents‑like and papers‑like runs.

**Example (patents‑like)**
```bash
./model_micro_macro 41 1.0 1.0 120000 0.0545 1.33 269781
```

**Example (papers‑like)**
```bash
./model_micro_macro 41 1.0 1.0 13000 0.175 0.89 295939
```

**Outputs (in `model_dir = data_simulations/`)**
- `model_T=..._nu=..._rho=..._N0=..._a=..._b=..._w0=....dat`
- `n_model_T=..._nu=..._rho=..._N0=..._a=..._b=..._w0=....dat`

### Reduce simulation files
**Notebook:** `module_3_reduce_data_sim.ipynb`

Reads the `.dat` files from `data_simulations/` and writes reduced CSVs into `output_dir`:
- Patents: `model_T=..._log.csv`, `model_T=..._lin.csv`, `model_T=..._freq.csv`
- Papers:  `outputs/model/model_T=..._log.csv`, `..._lin.csv`, `..._freq.csv`

## Module 4 — Comparison (Real Data vs Model)
**Notebook:** `module_4_comparison.ipynb`

**Inputs**
- Reduced real‑data series from `output_dir`
- Reduced simulation outputs from `data_simulations/` (via CSVs created above)

**What it does**
- Produces the 2×2 comparison panel (intrinsic time + real time).
- Compares frequency distributions (data vs simulation) for patents and papers.

## Module 5 — Branching Regimes (Model‑Only)

### Compile
```bash
gcc -o modelling_micro_macro_branching modelling_micro_macro_branching.c -lm
```

### Run (three regimes)
```bash
./modelling_micro_macro_branching 0 1000 1.0 1.0
./modelling_micro_macro_branching 1 600 1.0 1.0
./modelling_micro_macro_branching 2 44 1.0 1.0
```

**Outputs (in `model_dir = data_simulations/`)**
- `model_<label>_mode=<mode>_rho=..._nu=....dat`
- `n_model_<label>_mode=<mode>_rho=..._nu=....dat`

Then run `module_5_model.ipynb` to generate the figure comparing the three regimes.
