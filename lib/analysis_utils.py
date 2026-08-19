"""Derived datasets and model runners for the Supporting Information.

This module complements ``utils.py`` with two workflows:

1. the granularity datasets, which recount novelty at three resolutions
   (single elements, pairs, full combinations) and creators under both the
   focal-contributor and the all-contributor convention;
2. compilation, execution and scoring of the C model, used by the ``N0`` scan.

Event-level tables are written as Parquet next to a compact annual summary, so
that the SI notebooks can reload them without touching the raw data again.
"""

from __future__ import annotations

import ast
import json
import math
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from utils import (
    OPENALEX_SUBFIELDS,
    _paper_is_valid,
    _paper_year_mask,
    decimal_year,
    ensure_dir,
    normalize_keyword_combo,
    papers_mode_suffix,
    reduction_indices,
)


def parse_patent_inventors(raw_inventor_cell: object) -> list[str]:
    if isinstance(raw_inventor_cell, list):
        parsed = raw_inventor_cell
    else:
        parsed = None
        if isinstance(raw_inventor_cell, str):
            try:
                parsed = json.loads(raw_inventor_cell)
            except Exception:
                try:
                    parsed = ast.literal_eval(raw_inventor_cell)
                except Exception:
                    parsed = None
        if parsed is None:
            return []

    names: list[str] = []
    for entry in parsed:
        if isinstance(entry, str):
            try:
                entry = ast.literal_eval(entry)
            except Exception:
                continue
        if not isinstance(entry, dict):
            continue
        name = str(entry.get("name", "")).strip()
        if not name:
            continue
        if name.lower().startswith("the designation of the inventor"):
            continue
        names.append(name)
    return names


def parse_patent_ipc(raw_ipc_cell: object) -> list[str]:
    if isinstance(raw_ipc_cell, (list, tuple, np.ndarray)):
        values = list(raw_ipc_cell)
    elif isinstance(raw_ipc_cell, str):
        try:
            values = ast.literal_eval(raw_ipc_cell)
        except Exception:
            values = []
    else:
        values = []
    cleaned = sorted({str(v).strip() for v in values if str(v).strip()})
    return cleaned


def normalize_authors(authors: object) -> list[str]:
    if not isinstance(authors, (list, tuple, np.ndarray)):
        return []
    cleaned = []
    seen = set()
    for author in authors:
        name = str(author).strip()
        if not name or name in seen:
            continue
        seen.add(name)
        cleaned.append(name)
    return cleaned


def normalize_tokens(tokens: object) -> list[str]:
    if not isinstance(tokens, (list, tuple, np.ndarray)):
        return []
    return sorted({str(token).strip() for token in tokens if str(token).strip()})


def token_pairs(tokens: list[str]) -> list[tuple[str, str]]:
    if len(tokens) < 2:
        return []
    return list(combinations(tokens, 2))


def _relative_to_dir(value: object, base: Path) -> object:
    """Express ``value`` relative to ``base`` when it is a path underneath it.

    Used before persisting run inventories, so that a versioned CSV never
    records the absolute path of the machine that produced it.
    """
    if not isinstance(value, (str, Path)):
        return value
    try:
        return str(Path(value).resolve().relative_to(Path(base).resolve()))
    except ValueError:
        return str(value)


def _save_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))


def _metadata_path(output_path: Path) -> Path:
    return output_path.with_name(f"{output_path.stem}_meta.json")


@dataclass
class DatasetBundle:
    events_path: Path
    annual_path: Path
    metadata_path: Path


def _finalize_event_frame(payload: dict[str, list], output_path: Path) -> None:
    df = pd.DataFrame(payload)
    if "date" in df.columns:
        df["date"] = pd.to_datetime(df["date"])
    df.to_parquet(output_path, index=False)


def _finalize_annual_frame(rows: list[dict], output_path: Path) -> None:
    pd.DataFrame(rows).to_parquet(output_path, index=False)


def _annual_rows_from_counters(
    last_by_year: dict[int, dict[str, int | float]],
    active_first: dict[int, set[str]],
    active_all: dict[int, set[str]],
    team_size_sum: dict[int, int],
    event_count: dict[int, int],
) -> list[dict]:
    rows = []
    for year in sorted(last_by_year):
        base = dict(last_by_year[year])
        base["active_first_authors"] = len(active_first.get(year, set()))
        base["active_all_authors"] = len(active_all.get(year, set()))
        n_events = event_count.get(year, 0)
        base["mean_authors_per_event"] = (
            float(team_size_sum[year] / n_events) if n_events else np.nan
        )
        rows.append(base)
    return rows


def _load_papers_global_frame(
    papers_dir: str | Path,
    subfields: list[int],
    start_year: int,
    end_year: int,
) -> pd.DataFrame:
    papers_dir = Path(papers_dir)
    frames: list[pd.DataFrame] = []
    for idx, subfield in enumerate(subfields, start=1):
        path = papers_dir / f"merged_df_subfield_{subfield}.parquet"
        if not path.exists():
            continue
        print(f"[papers/global] [{idx}/{len(subfields)}] reading {path.name}")
        df = pd.read_parquet(path, columns=["date", "authors", "keywords", "topics"])
        mask_year = _paper_year_mask(df["date"], start_year, end_year)
        if not mask_year.any():
            continue
        df = df.loc[mask_year, ["date", "authors", "keywords", "topics"]].copy()
        df["date"] = pd.to_datetime(df["date"], errors="coerce")
        df = df.dropna(subset=["date"])
        frames.append(df)

    if not frames:
        return pd.DataFrame(columns=["date", "authors", "keywords", "topics"])

    print(f"[papers/global] concatenating {len(frames)} subfield chunks")
    full = pd.concat(frames, ignore_index=True)
    print(f"[papers/global] sorting {len(full):,} rows by date")
    full = full.sort_values("date", kind="stable").reset_index(drop=True)
    print(f"[papers/global] global frame ready: {len(full):,} rows")
    return full


def granularity_paths(
    output_dir: str | Path,
    filter_empty_keywords: bool = True,
) -> dict[str, DatasetBundle]:
    output_dir = ensure_dir(output_dir)
    suffix = papers_mode_suffix(filter_empty_keywords)
    return {
        "patents": DatasetBundle(
            events_path=output_dir / "patents_granularity.parquet",
            annual_path=output_dir / "patents_granularity_annual.parquet",
            metadata_path=output_dir / "patents_granularity_meta.json",
        ),
        "papers": DatasetBundle(
            events_path=output_dir / f"papers_granularity_{suffix}.parquet",
            annual_path=output_dir / f"papers_granularity_{suffix}_annual.parquet",
            metadata_path=output_dir / f"papers_granularity_{suffix}_meta.json",
        ),
    }


def _append_granularity_payload(
    payload: dict[str, list],
    ts: pd.Timestamp,
    t: int,
    d_single: int,
    d_pair: int,
    d_combo: int,
    dw_first: int,
    dw_all: int,
    n_authors_event: int,
    n_tokens_event: int,
) -> None:
    payload["date"].append(ts)
    payload["date_ord"].append(int(ts.toordinal()))
    payload["tau"].append(float(decimal_year(ts)))
    payload["year"].append(int(ts.year))
    payload["t"].append(t)
    payload["D_single"].append(d_single)
    payload["D_pair"].append(d_pair)
    payload["D_combo"].append(d_combo)
    payload["Dw_first"].append(dw_first)
    payload["Dw_all"].append(dw_all)
    payload["n_authors_event"].append(n_authors_event)
    payload["n_tokens_event"].append(n_tokens_event)


def build_patents_granularity_dataset(
    patents_aug_dir: str | Path,
    output_dir: str | Path,
    start_year: int = 1980,
    end_year: int = 2020,
) -> DatasetBundle:
    bundle = granularity_paths(output_dir)["patents"]

    seen_single: set[str] = set()
    seen_pair: set[tuple[str, str]] = set()
    seen_combo: set[tuple[str, ...]] = set()
    seen_first: set[str] = set()
    seen_all: set[str] = set()

    active_first: dict[int, set[str]] = defaultdict(set)
    active_all: dict[int, set[str]] = defaultdict(set)
    team_size_sum: dict[int, int] = defaultdict(int)
    event_count: dict[int, int] = defaultdict(int)
    last_by_year: dict[int, dict[str, int | float]] = {}

    payload = {
        "date": [],
        "date_ord": [],
        "tau": [],
        "year": [],
        "t": [],
        "D_single": [],
        "D_pair": [],
        "D_combo": [],
        "Dw_first": [],
        "Dw_all": [],
        "n_authors_event": [],
        "n_tokens_event": [],
    }

    t = d_single = d_pair = d_combo = dw_first = dw_all = 0

    for path in sorted(Path(patents_aug_dir).glob("*.csv")):
        if not path.stem.isdigit():
            continue
        year = int(path.stem)
        if year < start_year or year > end_year:
            continue
        print(f"[granularity/patents] reading {path.name}")
        df = pd.read_csv(path, usecols=["publication.date", "inventor", "ipc"])
        df["publication.date"] = pd.to_datetime(df["publication.date"], errors="coerce")
        df = df.dropna(subset=["publication.date"]).sort_values("publication.date", kind="stable")

        for row in df.itertuples(index=False):
            authors = parse_patent_inventors(row.inventor)
            if not authors:
                continue
            tokens = parse_patent_ipc(row.ipc)
            combo = tuple(tokens)
            pairs = token_pairs(tokens)
            first_author = authors[0]

            t += 1
            for token in tokens:
                if token not in seen_single:
                    seen_single.add(token)
                    d_single += 1
            for pair in pairs:
                if pair not in seen_pair:
                    seen_pair.add(pair)
                    d_pair += 1
            if combo not in seen_combo:
                seen_combo.add(combo)
                d_combo += 1
            if first_author not in seen_first:
                seen_first.add(first_author)
                dw_first += 1
            new_all = 0
            for author in authors:
                if author not in seen_all:
                    seen_all.add(author)
                    new_all += 1
            dw_all += new_all

            ts = pd.Timestamp(row[0])
            year_val = int(ts.year)
            active_first[year_val].add(first_author)
            active_all[year_val].update(authors)
            team_size_sum[year_val] += len(authors)
            event_count[year_val] += 1

            _append_granularity_payload(
                payload,
                ts,
                t,
                d_single,
                d_pair,
                d_combo,
                dw_first,
                dw_all,
                len(authors),
                len(tokens),
            )
            last_by_year[year_val] = {
                "year": year_val,
                "tau": year_val,
                "delta_t": event_count[year_val],
                "t": t,
                "D_single": d_single,
                "D_pair": d_pair,
                "D_combo": d_combo,
                "Dw_first": dw_first,
                "Dw_all": dw_all,
            }

    _finalize_event_frame(payload, bundle.events_path)
    _finalize_annual_frame(
        _annual_rows_from_counters(last_by_year, active_first, active_all, team_size_sum, event_count),
        bundle.annual_path,
    )
    _save_json(
        bundle.metadata_path,
        {
            "source": str(patents_aug_dir),
            "n_events": int(t),
            "final_D_single": int(d_single),
            "final_D_pair": int(d_pair),
            "final_D_combo": int(d_combo),
            "final_Dw_first": int(dw_first),
            "final_Dw_all": int(dw_all),
        },
    )
    return bundle


def build_papers_granularity_dataset(
    papers_dir: str | Path,
    output_dir: str | Path,
    start_year: int = 1980,
    end_year: int = 2020,
    subfields: list[int] | None = None,
    filter_empty_keywords: bool = True,
) -> DatasetBundle:
    bundle = granularity_paths(output_dir, filter_empty_keywords=filter_empty_keywords)["papers"]

    if subfields is None:
        subfields = OPENALEX_SUBFIELDS

    seen_single: set[str] = set()
    seen_pair: set[tuple[str, str]] = set()
    seen_combo: set[str] = set()
    seen_first: set[str] = set()
    seen_all: set[str] = set()

    active_first: dict[int, set[str]] = defaultdict(set)
    active_all: dict[int, set[str]] = defaultdict(set)
    team_size_sum: dict[int, int] = defaultdict(int)
    event_count: dict[int, int] = defaultdict(int)
    last_by_year: dict[int, dict[str, int | float]] = {}

    payload = {
        "date": [],
        "date_ord": [],
        "tau": [],
        "year": [],
        "t": [],
        "D_single": [],
        "D_pair": [],
        "D_combo": [],
        "Dw_first": [],
        "Dw_all": [],
        "n_authors_event": [],
        "n_tokens_event": [],
    }

    t = d_single = d_pair = d_combo = dw_first = dw_all = 0
    progress_step = 250_000
    df = _load_papers_global_frame(
        papers_dir=papers_dir,
        subfields=subfields,
        start_year=start_year,
        end_year=end_year,
    )
    print(f"[granularity/papers] processing globally ordered rows: {len(df):,}")

    for idx, (date_value, authors_raw, keywords_raw, topics_raw) in enumerate(zip(
        df["date"],
        df["authors"],
        df["keywords"],
        df["topics"],
    ), start=1):
        if idx == 1 or idx % progress_step == 0 or idx == len(df):
            print(
                f"[granularity/papers] progress {idx:,}/{len(df):,} | "
                f"kept t={t:,} D_single={d_single:,} D_pair={d_pair:,} "
                f"D_combo={d_combo:,} Dw_first={dw_first:,} Dw_all={dw_all:,}"
            )
        if not _paper_is_valid(
            date_value,
            authors_raw,
            keywords_raw,
            topics_raw,
            require_keywords=filter_empty_keywords,
        ):
            continue
        authors = normalize_authors(authors_raw)
        if not authors:
            continue
        tokens = normalize_tokens(keywords_raw)
        combo = normalize_keyword_combo(tokens) if tokens else ""
        pairs = token_pairs(tokens)

        ts = pd.to_datetime(date_value, errors="coerce")
        if pd.isna(ts):
            continue

        t += 1
        for token in tokens:
            if token not in seen_single:
                seen_single.add(token)
                d_single += 1
        for pair in pairs:
            if pair not in seen_pair:
                seen_pair.add(pair)
                d_pair += 1
        if combo and combo not in seen_combo:
            seen_combo.add(combo)
            d_combo += 1
        first_author = authors[0]
        if first_author not in seen_first:
            seen_first.add(first_author)
            dw_first += 1
        new_all = 0
        for author in authors:
            if author not in seen_all:
                seen_all.add(author)
                new_all += 1
        dw_all += new_all

        year_val = int(ts.year)
        active_first[year_val].add(first_author)
        active_all[year_val].update(authors)
        team_size_sum[year_val] += len(authors)
        event_count[year_val] += 1

        _append_granularity_payload(
            payload,
            ts,
            t,
            d_single,
            d_pair,
            d_combo,
            dw_first,
            dw_all,
            len(authors),
            len(tokens),
        )
        last_by_year[year_val] = {
            "year": year_val,
            "tau": year_val,
            "delta_t": event_count[year_val],
            "t": t,
            "D_single": d_single,
            "D_pair": d_pair,
            "D_combo": d_combo,
            "Dw_first": dw_first,
            "Dw_all": dw_all,
        }

    _finalize_event_frame(payload, bundle.events_path)
    _finalize_annual_frame(
        _annual_rows_from_counters(last_by_year, active_first, active_all, team_size_sum, event_count),
        bundle.annual_path,
    )
    _save_json(
        bundle.metadata_path,
        {
            "source": str(papers_dir),
            "subfields": sorted(subfields),
            "filter_empty_keywords": bool(filter_empty_keywords),
            "n_events": int(t),
            "final_D_single": int(d_single),
            "final_D_pair": int(d_pair),
            "final_D_combo": int(d_combo),
            "final_Dw_first": int(dw_first),
            "final_Dw_all": int(dw_all),
        },
    )
    return bundle


def build_granularity_datasets(
    output_dir: str | Path,
    patents_aug_dir: str | Path,
    papers_dir: str | Path,
    start_year: int = 1980,
    end_year: int = 2020,
    subfields: list[int] | None = None,
    filter_empty_keywords: bool = True,
) -> dict[str, DatasetBundle]:
    return {
        "patents": build_patents_granularity_dataset(
            patents_aug_dir=patents_aug_dir,
            output_dir=output_dir,
            start_year=start_year,
            end_year=end_year,
        ),
        "papers": build_papers_granularity_dataset(
            papers_dir=papers_dir,
            output_dir=output_dir,
            start_year=start_year,
            end_year=end_year,
            subfields=subfields,
            filter_empty_keywords=filter_empty_keywords,
        ),
    }


def load_granularity_datasets(
    output_dir: str | Path,
    filter_empty_keywords: bool = True,
) -> dict[str, dict[str, pd.DataFrame]]:
    bundles = granularity_paths(output_dir, filter_empty_keywords=filter_empty_keywords)
    out = {}
    for dataset, bundle in bundles.items():
        out[dataset] = {
            "events": pd.read_parquet(bundle.events_path),
            "annual": pd.read_parquet(bundle.annual_path),
        }
    return out


def compile_model_binary(
    source_dir: str | Path,
    binary_name: str = "model_micro_macro",
    source_name: str = "model_micro_macro.c",
    extra_flags: list[str] | None = None,
) -> Path:
    source_dir = Path(source_dir)
    binary_path = source_dir / binary_name
    # gnu11, not c11: srand48/lrand48 are POSIX, and under strict c11 they are
    # only implicitly declared. The output is identical either way, but the
    # implicit declarations are a portability hazard on other toolchains.
    cmd = ["gcc", "-O3", "-std=gnu11", source_name, "-lm", "-o", binary_name]
    if extra_flags:
        cmd.extend(extra_flags)
    print(f"[model] compiling: {' '.join(cmd)} (cwd={source_dir})")
    subprocess.run(cmd, cwd=source_dir, check=True)
    return binary_path


def expected_model_output_paths(
    workdir: str | Path,
    T: int,
    nu: float,
    rho: float,
    N0: int,
    a: float,
    b: float,
    w0: float,
    output_stem: str = "model",
) -> dict[str, Path]:
    workdir = Path(workdir)
    prefix = f"{output_stem}_T={T}_nu={nu:.6g}_rho={rho:.6g}_N0={N0}_a={a:.6g}_b={b:.6g}_w0={w0:.6g}"
    return {
        "trajectory": workdir / "data_simulations" / f"{prefix}.dat",
        "frequency": workdir / "data_simulations" / f"n_{prefix}.dat",
    }


def run_model_simulation(
    binary_path: str | Path,
    workdir: str | Path,
    T: int,
    nu: float,
    rho: float,
    N0: int,
    a: float,
    b: float,
    w0: float,
    output_stem: str = "model",
    verbose: bool = True,
) -> dict[str, object]:
    binary_path = Path(binary_path)
    workdir = ensure_dir(workdir)
    expected = expected_model_output_paths(
        workdir,
        T,
        nu,
        rho,
        N0,
        a,
        b,
        w0,
        output_stem=output_stem,
    )

    ensure_dir(workdir / "data_simulations")
    cmd = [
        str(binary_path),
        str(T),
        str(nu),
        str(rho),
        str(N0),
        str(a),
        str(b),
        str(w0),
    ]
    print(f"[model] running: {' '.join(cmd)} (cwd={workdir})")
    run_kwargs = {"cwd": workdir, "check": True}
    if not verbose:
        run_kwargs["stdout"] = subprocess.DEVNULL
        run_kwargs["stderr"] = subprocess.DEVNULL
    subprocess.run(cmd, **run_kwargs)
    return {"N0": N0, **expected}


def reduce_simulation_run(
    dat_path: str | Path,
    output_dir: str | Path,
    model_prefix: str,
    n_points: int = 10_000,
) -> dict[str, int]:
    """Reduce one raw model trajectory to the three CSVs Module 4 reads.

    The raw ``.dat`` file has one row per extraction and columns
    ``tau, t, D, estr, num_estr, w``; it is far too large to plot directly.
    Writes into ``output_dir``:

    - ``<model_prefix>_log.csv`` and ``<model_prefix>_lin.csv`` — the
      trajectory subsampled on the log and linear grids;
    - ``<model_prefix>_freq.csv`` — the rank-frequency table of the extracted
      elements.

    Returns the end-of-run totals, for the printed sanity check.
    """

    dat_path = Path(dat_path)
    output_dir = ensure_dir(output_dir)

    with open(dat_path, "r") as handle:
        lines = handle.readlines()
    try:
        data = np.loadtxt(lines)
    except ValueError as exc:
        # A run interrupted mid-write leaves a truncated final row.
        print(f"Warning '{exc}', dropping last corrupted line.")
        data = np.loadtxt(lines[:-1])

    tau, t, D, estr, num_estr, w = (data[:, i] for i in range(6))

    freq = np.array(sorted(Counter(estr).values(), reverse=True))
    pd.DataFrame({
        "rank": np.arange(1, len(freq) + 1),
        "frequency": freq,
    }).to_csv(output_dir / f"{model_prefix}_freq.csv", index=False)

    log_indices, lin_indices = reduction_indices(len(tau), n_points=n_points)
    for suffix, idx in (("log", log_indices), ("lin", lin_indices)):
        pd.DataFrame({
            "tau": tau[idx],
            "t": t[idx],
            "D": D[idx],
            "num_estr": num_estr[idx],
            "w": w[idx],
        }).to_csv(output_dir / f"{model_prefix}_{suffix}.csv", index=False)

    return {
        "rows": len(tau),
        "events": int(t[-1]),
        "explorers": int(w.max()),
        "novelties": int(D.max()),
    }


def load_model_trajectory(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=["tau", "t", "D", "estr", "num_estr", "Dw"],
    )
    return df


def model_annual_summary(
    model_df: pd.DataFrame,
    year0: int = 1980,
) -> pd.DataFrame:
    if model_df.empty:
        return pd.DataFrame(columns=["year", "tau", "delta_t", "t", "D", "Dw"])
    tmp = model_df.copy()
    tmp["tau_year"] = np.floor(tmp["tau"]).astype(int)
    grouped = tmp.groupby("tau_year", sort=True)
    out = grouped.agg(
        t=("t", "last"),
        D=("D", "last"),
        Dw=("Dw", "last"),
    ).reset_index()
    out["year"] = year0 + out["tau_year"]
    out["tau"] = out["year"]
    out["delta_t"] = out["t"].diff().fillna(out["t"]).astype(float)
    return out[["year", "tau", "delta_t", "t", "D", "Dw"]]


def score_model_vs_empirical(
    model_annual: pd.DataFrame,
    empirical_annual: pd.DataFrame,
    columns: tuple[str, ...] = ("t", "D", "Dw"),
    combined_columns: tuple[str, ...] | None = None,
    start_year: int | None = None,
) -> dict[str, float]:
    left = empirical_annual.copy()
    right = model_annual.copy()
    if start_year is not None:
        left = left[left["year"] >= start_year]
        right = right[right["year"] >= start_year]
    merged = left.merge(right, on="year", suffixes=("_emp", "_mod"))
    if merged.empty:
        return {f"rmse_{col}": np.nan for col in columns} | {"rmse_log_total": np.nan}

    if combined_columns is None:
        combined_columns = columns
    unknown_combined = set(combined_columns) - set(columns)
    if unknown_combined:
        raise ValueError(
            f"combined_columns must be included in columns; unknown: {unknown_combined}"
        )

    metrics: dict[str, float] = {}
    total_log = 0.0
    used = 0
    for col in columns:
        emp = merged[f"{col}_emp"].to_numpy(dtype=float)
        mod = merged[f"{col}_mod"].to_numpy(dtype=float)
        rmse = float(np.sqrt(np.mean((mod - emp) ** 2)))
        metrics[f"rmse_{col}"] = rmse
        if col in combined_columns and np.all(emp > 0) and np.all(mod > 0):
            total_log += float(np.mean((np.log(mod) - np.log(emp)) ** 2))
            used += 1
    metrics["rmse_log_total"] = float(math.sqrt(total_log / used)) if used else np.nan
    return metrics


def score_intrinsic_model_vs_empirical(
    model_df: pd.DataFrame,
    empirical_df: pd.DataFrame,
    columns: tuple[str, ...] = ("D", "Dw"),
    n_points: int = 5000,
) -> dict[str, float]:
    """Score intrinsic-time trajectories on a dense, uniformly spaced t grid."""
    t_min = max(
        1.0,
        float(model_df.loc[model_df["t"] > 0, "t"].min()),
        float(empirical_df.loc[empirical_df["t"] > 0, "t"].min()),
    )
    t_max = min(float(model_df["t"].max()), float(empirical_df["t"].max()))
    if t_max <= t_min:
        return {f"rmse_{col}": np.nan for col in columns} | {"rmse_log_total": np.nan}

    sample_t = np.unique(
        np.rint(np.linspace(t_min, t_max, num=n_points)).astype(np.int64)
    ).astype(float)
    metrics: dict[str, float] = {}
    total_log = 0.0
    used = 0
    for col in columns:
        empirical_values = np.interp(
            sample_t,
            empirical_df["t"].to_numpy(dtype=float),
            empirical_df[col].to_numpy(dtype=float),
        )
        model_values = np.interp(
            sample_t,
            model_df["t"].to_numpy(dtype=float),
            model_df[col].to_numpy(dtype=float),
        )
        metrics[f"rmse_{col}"] = float(
            np.sqrt(np.mean((model_values - empirical_values) ** 2))
        )
        positive = (empirical_values > 0) & (model_values > 0)
        if np.any(positive):
            total_log += float(
                np.mean(
                    (
                        np.log(model_values[positive])
                        - np.log(empirical_values[positive])
                    )
                    ** 2
                )
            )
            used += 1
    metrics["rmse_log_total"] = float(math.sqrt(total_log / used)) if used else np.nan
    return metrics


def run_and_score_model_grid(
    binary_path: str | Path,
    workdir: str | Path,
    base_params: dict[str, float | int],
    N0_values: Iterable[int],
    empirical_annual: pd.DataFrame,
    empirical_intrinsic: pd.DataFrame | None = None,
    intrinsic_points: int = 5000,
    start_year: int | None = None,
    combined_columns: tuple[str, ...] | None = None,
    output_stem: str = "model",
    cleanup_outputs: bool = True,
    save_results_csv: str | Path | None = None,
    verbose: bool = True,
) -> pd.DataFrame:
    requested_n0_values = [int(x) for x in N0_values]
    if save_results_csv is not None:
        save_results_csv = Path(save_results_csv)

    rows: list[dict[str, object]] = []
    year0 = int(empirical_annual["year"].min())
    for N0 in requested_n0_values:
        N0 = int(N0)
        params = dict(base_params)
        params["N0"] = N0
        out = run_model_simulation(
            binary_path=binary_path,
            workdir=workdir,
            T=int(params["T"]),
            nu=float(params["nu"]),
            rho=float(params["rho"]),
            N0=int(params["N0"]),
            a=float(params["a"]),
            b=float(params["b"]),
            w0=float(params["w0"]),
            output_stem=output_stem,
            verbose=verbose,
        )
        model_df = load_model_trajectory(out["trajectory"])
        if empirical_intrinsic is not None:
            metrics = score_intrinsic_model_vs_empirical(
                model_df,
                empirical_intrinsic,
                columns=combined_columns or ("D", "Dw"),
                n_points=intrinsic_points,
            )
        else:
            annual_df = model_annual_summary(model_df, year0=year0)
            metrics = score_model_vs_empirical(
                annual_df,
                empirical_annual,
                combined_columns=combined_columns,
                start_year=start_year,
            )
        row = dict(params)
        row["trajectory_path"] = str(out["trajectory"])
        row["frequency_path"] = str(out["frequency"])
        row |= metrics
        rows.append(row)
        if cleanup_outputs:
            for key in ("trajectory", "frequency"):
                path = Path(out[key])
                if path.exists():
                    path.unlink()
        print(
            "[model] finished N0="
            f"{int(params['N0'])} | combined log-RMSE={row['rmse_log_total']:.6f}"
        )

    results = pd.DataFrame(rows)
    if not results.empty:
        results = results.sort_values("N0").drop_duplicates(subset=["N0"], keep="last")
    if save_results_csv is not None:
        ensure_dir(save_results_csv.parent)
        save_df = results.copy()
        # Persist the trajectory/frequency locations relative to the results
        # file: the grid CSV is versioned, so it must not carry the absolute
        # path of whoever ran the scan. The columns are informational; nothing
        # reads them back (the .dat files are removed by cleanup_outputs).
        for _col in ("trajectory_path", "frequency_path"):
            if _col in save_df.columns:
                save_df[_col] = [
                    _relative_to_dir(value, save_results_csv.parent)
                    for value in save_df[_col]
                ]
        save_df.to_csv(save_results_csv, index=False)
        print(f"[model] saved grid results to {save_results_csv}")
    return results
