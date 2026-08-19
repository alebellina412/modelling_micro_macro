# The gold-rush effect: how innovation speeds up

Code and analysis pipeline for the paper. Everything the manuscript reports is
produced here: ten of its eleven figures (all but the schematic), both tables,
and the numbers quoted in the text and in the Supporting Information.

The pipeline runs in four stages:

1. **build** the empirical event series from the raw patent and paper archives;
2. **calibrate** the model against those series and pick its one free
   parameter, the initial knowledge stock `N0`;
3. **simulate** the model and reduce the output;
4. **compare** data and model, and produce the figures and tables.

## Calibration

Three of the model parameters — `a`, `b` and `w0` — are fixed by the empirical
fits of Module 2. The fourth, the initial knowledge stock `N0`, is the one the
fits do not pin down, and the scan of Module 4 selects it. The values used
throughout (`T = 41`, `nu = rho = 1` in both cases):

| | `a` | `b` | `w0` | `N0` |
|---|---|---|---|---|
| patents EPO | 0.0545 | 1.33 | 269,781 | 100,500 |
| papers OA | 0.167 | 0.71 | 212,313 | 49,000 |

The three fitted values are carried over by hand: Module 2 prints them, and they
appear as literals in the parameter cells downstream.

## Requirements

Python 3.10+ and a C compiler. Install the Python side with:

```bash
pip install -r requirements.txt
```

The figure style is fixed centrally in `lib/plot_style.py`; notebooks do not
override it locally.

## Getting the data

The raw archives are third-party data and are not redistributed here. Create a
`data/` directory at the repository root and place them as follows.

- **`data/patents_EPO/`** — one CSV per year, named `YYYY.csv` (`1980.csv` …
  `2020.csv`), each with the columns `publication.date` (`YYYYMMDD`),
  `inventor` (a list-like string of dicts, each with a `name`) and `ipc` (a
  list-like string of IPC codes).
- **`data/papers_openalex/`** — one parquet per OpenAlex subfield, named
  `merged_df_subfield_{ID}.parquet`, each with the columns `date`, `topics`,
  `keywords` and `authors`. The eighteen subfield IDs are listed in the
  parameter cell of `module_1_papers.ipynb`.

If the archives already live elsewhere on your machine, symlink them instead of
copying:

```bash
mkdir -p data
ln -s /path/to/patents_EPO      data/patents_EPO
ln -s /path/to/papers_openalex  data/papers_openalex
```

`data/` is excluded from version control, as are the two directories the
pipeline generates: `outputs/` (derived data) and `data_simulations/` (raw
simulation trajectories). Only figures, tables and the small run inventories in
`sim_runs/` are versioned.

## How to run it

Run the notebooks **from the repository root**, in this order. Each one starts
with a bootstrap cell that walks up to the directory containing `lib/` and puts
it on `sys.path`, so nothing is imported from outside the repository; the C
programs write relative to the working directory, so they too must be launched
from the root.

The module numbers are the run order.

| # | step | needs | produces |
|---|---|---|---|
| 1a | `module_1_patents.ipynb` | raw patents | patent series in `outputs/` |
| 1b | `module_1_papers.ipynb` | raw papers | paper series in `outputs/` |
| 2 | `module_2_real_data.ipynb` | 1a, 1b | the fits `a`, `b`, `w0`; `fig_1`, `fig_6`, `fig_SI1` |
| 3 | `module_3_trajectories.ipynb` | raw papers, 1a | `outputs/trajectories/` |
| 4 | `module_4_n0_search.ipynb` | 3 | `N0` scan and selected run, for both datasets |
| 5 | `module_5_reduce_data_sim.ipynb` | 4 | reduced model series |
| 6 | `module_6_comparison.ipynb` | 1a, 1b, 5 | `fig_3`, `fig_4` |
| 7 | `modelling_micro_macro_branching` + `module_7_model.ipynb` | — | `fig_5` |
| 8 | `si_build_granularity.ipynb` | raw papers, 1a | `outputs/granularity/` |
| 9 | `si_model_selection.ipynb` | 3 | `fig_SI3`, both tables |
| 10 | `si_granularity.ipynb` | 8 | `fig_SI4` |
| 11 | `si_contributors.ipynb` | 8 | `fig_SI2` |
| 12 | `si_n0_sensitivity.ipynb` | 4 | `fig_SI5` |

Module 2 comes before the calibration because it is where `a`, `b` and `w0` are
fitted; Module 7 depends on nothing and can be run at any point.

Steps 1, 3 and 8 are the long ones: they stream the full archives. Step 4 runs
the model 54 times, a few seconds each, but re-reads every trajectory to score
it. Everything else works from cached intermediates and takes about three
minutes in total on a desktop machine, `module_2_real_data.ipynb` being the
slowest at around 80 seconds.

Step 7 needs the three branching runs first, from the repository root:

```bash
gcc -O3 -std=gnu11 -o modelling_micro_macro_branching modelling_micro_macro_branching.c -lm
./modelling_micro_macro_branching 0 1000 1.0 1.0   # logarithmic
./modelling_micro_macro_branching 1  600 1.0 1.0   # simple (linear) branching
./modelling_micro_macro_branching 2   44 1.0 1.0   # singularity
```

### Reproducibility

The C programs seed with `srand48(1)`, so the simulations are deterministic:
the same parameters give the same trajectory, byte for byte, and the figures
come out byte-identical to the ones in the manuscript.

Compile with `-std=gnu11`, not `-std=c11`: `srand48` and `lrand48` are POSIX,
and under strict C11 they are only implicitly declared. The output is the same
either way on glibc/x86-64, but the implicit declarations are a hazard on other
toolchains.

## Figures and tables

`figures/` and `tables/` hold exactly what the manuscript includes, under the
same names. `fig_2.pdf` is the schematic of the model and is the only figure
not produced by this repository.

| file | notebook | where it appears |
|---|---|---|
| `fig_1.pdf` | `module_2_real_data.ipynb` | main, Fig. 1 |
| `fig_3.pdf` | `module_6_comparison.ipynb` | main, Fig. 3 |
| `fig_4.pdf` | `module_6_comparison.ipynb` | main, Fig. 4 |
| `fig_5.pdf` | `module_7_model.ipynb` | main, Fig. 5 |
| `fig_6.pdf` | `module_2_real_data.ipynb` | main, Fig. 6 |
| `fig_SI1.pdf` | `module_2_real_data.ipynb` | SI, Fig. S1 |
| `fig_SI2.pdf` | `si_contributors.ipynb` | SI, Fig. S2 |
| `fig_SI3.pdf` | `si_model_selection.ipynb` | SI, Fig. S3 |
| `fig_SI4.pdf` | `si_granularity.ipynb` | SI, Fig. S4 |
| `fig_SI5.pdf` | `si_n0_sensitivity.ipynb` | SI, Fig. S5 |
| `table_intrinsic_fits.tex` | `si_model_selection.ipynb` | SI, Table S1 |
| `table_branching_fits.tex` | `si_model_selection.ipynb` | SI, Table S2 |

## Layout

```
lib/                  shared Python modules, imported by every notebook
module_*.ipynb        the main pipeline
si_*.ipynb            the Supporting Information analyses
*.c, polya_adj.h      the model
figures/, tables/     manuscript figures and tables (versioned)
sim_runs/             model-run inventories (versioned; trajectories are not)
data/                 raw inputs, you provide these       (not versioned)
outputs/              derived data                        (not versioned)
data_simulations/     raw simulation trajectories         (not versioned)
```

### `lib/`

- `plot_style.py` — the frozen figure style `pnas-main-v1`, shared by the main
  and the SI figures;
- `utils.py` — trajectory building, intrinsic-time and branching fits, binned
  effective exponents;
- `analysis_utils.py` — the granularity and contributor datasets, and the
  compilation, execution and scoring of the C model;
- `module1_optimized.py` — the Module 1 papers builder (the patents
  builder lives inline in `module_1_patents.ipynb`).

### The modules

**Module 1 — real-data preparation.** `module_1_patents.ipynb` and
`module_1_papers.ipynb` turn the raw archives into the event series. For
patents, novelty is the set of IPC codes of a filing and the explorer is its
first inventor; for papers, novelty is the keyword combination and the explorer
is the first author.

An event only counts if it can carry a novelty. A publication with no keyword
cannot form a combination, so `module_1_papers.ipynb` runs with
`filter_empty_keywords = True` and leaves those papers out of the sequence; the
patent archive needs no equivalent filter, since every filing in it carries at
least one IPC code.

Both notebooks drop first authors with implausibly many events in a single year
(institutional filers rather than individuals), then accumulate `t`, `D` and
`D_w` over the event sequence and write the reduced series, the frequency tables
and the per-author event timelines into `outputs/`. The patents notebook carries
its own code; the papers one delegates to `lib/module1_optimized.py`, a
streaming builder that stages through SQLite instead of holding the archive in
memory.

**Module 2 — fits and rates.** Fits `D_w = b D` and `dt/dtau` against `w` to
obtain `b`, `a` and `w0`, and produces the empirical figures of the main text
together with the inter-event time and career-length statistics of the SI.

**Module 3 — event trajectories.** One row per event, with the cumulative `t`,
`D` and `D_w`, for both datasets. This is the table the model is calibrated
against.

**Module 4 — the initial knowledge stock.** `N0` is the one parameter the fits
do not pin down. The notebook compiles the model, runs it over a grid of `N0`
values for both datasets, scores each run against the empirical trajectory by
the combined log-RMSE on `D` and `D_w`, and re-runs the best one. Each
trajectory is deleted right after it is scored; only the inventory
(`n0_grid_results.csv`, `best_n0_run.json`) and the winning run are kept, under
`sim_runs/<dataset>/`.

**Module 5 — reducing the simulation output.** The model writes one row per
extraction, tens of millions of rows per run. This notebook subsamples each run
onto a log and a linear grid and extracts the rank-frequency table.

**Module 6 — comparison.** Data against model, in natural and intrinsic time,
and the novelty frequency distributions.

**Module 7 — branching regimes.** The model alone, under the three branching
laws: logarithmic, linear, and the singular case.

### The Supporting Information notebooks

Two of the four read a cache that has to be built first. `si_granularity.ipynb`
and `si_contributors.ipynb` need `outputs/granularity/`, written by
`si_build_granularity.ipynb`: the same events as Module 3, but with novelty
counted at three resolutions and creators counted under both contributor
conventions. `si_model_selection.ipynb` reads the Module 3 trajectory cache
instead, and `si_n0_sensitivity.ipynb` needs only `sim_runs/`.

- `si_model_selection.ipynb` — compares the candidate intrinsic-time laws for
  `D(t)` (linear, `t/log t`, power law) and the candidate branching laws for
  `D_w(D)` (linear, logarithmic, power law, exponential), selecting on linear
  RMSE, and writes the two fit tables.
- `si_granularity.ipynb` — the same two laws at three novelty resolutions:
  single elements, pairs, full combinations.
- `si_contributors.ipynb` — the cumulative all-contributor count against the
  cumulative focal-contributor count. The slope of the linear fit is the
  effective rescaling factor between the two conventions; it is not the mean
  team size, because contributors recur across events.
- `si_n0_sensitivity.ipynb` — the combined log-RMSE as a function of `N0`,
  showing that the selected value sits in a flat basin rather than at a sharp
  minimum.

## The model

```bash
gcc -O3 -std=gnu11 -o model_micro_macro model_micro_macro.c -lm

./model_micro_macro T nu rho N0 a b w0
```

Use `T = 41` for both patents-like and papers-like runs:

```bash
./model_micro_macro 41 1.0 1.0 100500 0.0545 1.33 269781   # patents-like
./model_micro_macro 41 1.0 1.0  49000 0.167  0.71 212313   # papers-like
```

Both write a trajectory and a frequency file into `data_simulations/`.
`module_4_n0_search.ipynb` compiles and drives the model itself, so for the
published pipeline these commands are only needed to reproduce a single run by
hand.

## License

Released under the MIT License; see [`LICENSE`](LICENSE). The raw archives under
`data/` are not covered by it: they are third-party data, are not redistributed
here, and keep the terms of their own providers.
