"""Trajectory building, scaling diagnostics and fits.

This module works from full trajectory files saved on disk. The trajectory
files are ordinary Parquet datasets with one row per event and cumulative
columns already materialized:

- patents_trajectory.parquet
- papers_trajectory.parquet

For papers, the builder uses a disk-backed SQLite staging database so that the
final trajectory is exact but does not require keeping the full rich-text raw
data in memory. The switch is `filter_empty_keywords`:

- True: count only papers with valid keyword combinations
- False: reproduce the public github pipeline logic, where papers without
  keywords still contribute to t and Dw but not to D
"""

from __future__ import annotations

import json
import math
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.stats as stats
import statsmodels.api as sm

OPENALEX_SUBFIELDS = [14, 25, 37, 39, 54, 110, 114, 131, 142, 164, 172, 191, 193, 206, 210, 216, 246, 248]


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def reduction_indices(n: int, n_points: int = 10_000) -> tuple[np.ndarray, np.ndarray]:
    """Log- and linearly-spaced sample indices for a series of length ``n``.

    Every long series in the pipeline — empirical and simulated alike — is
    stored twice, subsampled this way: the log grid resolves the early growth,
    the linear grid the late one.
    """

    n_points = min(n_points, n)
    log_indices = np.unique(np.logspace(0, np.log10(n), n_points, dtype=int) - 1)
    lin_indices = np.linspace(0, n - 1, n_points, dtype=int)
    return log_indices.astype(np.int64), lin_indices.astype(np.int64)


def decimal_year(ts: pd.Timestamp) -> float:
    day_of_year = ts.timetuple().tm_yday
    return float(ts.year + (day_of_year - 1) / 365.25)


def normalize_keyword_combo(keywords: Iterable[object]) -> str:
    normalized = sorted({str(item).strip() for item in keywords if str(item).strip()})
    return "\x1f".join(normalized)


def papers_mode_suffix(filter_empty_keywords: bool) -> str:
    return "filter_empty_keywords_true" if filter_empty_keywords else "filter_empty_keywords_false"


def trajectory_paths(
    trajectory_dir: str | Path,
    filter_empty_keywords: bool = True,
) -> dict[str, Path]:
    trajectory_dir = ensure_dir(trajectory_dir)
    suffix = papers_mode_suffix(filter_empty_keywords)
    return {
        "patents": trajectory_dir / "patents_trajectory.parquet",
        "papers": trajectory_dir / f"papers_trajectory_{suffix}.parquet",
    }


def metadata_path(trajectory_dir: str | Path, prefix: str) -> Path:
    return ensure_dir(trajectory_dir) / f"{prefix}_trajectory_meta.json"


def papers_metadata_path(trajectory_dir: str | Path, filter_empty_keywords: bool) -> Path:
    suffix = papers_mode_suffix(filter_empty_keywords)
    return ensure_dir(trajectory_dir) / f"papers_trajectory_{suffix}_meta.json"


@dataclass
class Trajectory:
    name: str
    df: pd.DataFrame

    @property
    def n(self) -> int:
        return int(len(self.df))

    @property
    def t(self) -> np.ndarray:
        return self.df["t"].to_numpy()

    @property
    def D(self) -> np.ndarray:
        return self.df["D"].to_numpy()

    @property
    def Dw(self) -> np.ndarray:
        return self.df["Dw"].to_numpy()

    @property
    def tau(self) -> np.ndarray:
        return self.df["tau"].to_numpy()

    @property
    def year(self) -> np.ndarray:
        return self.df["year"].to_numpy()

    @property
    def date_ord(self) -> np.ndarray:
        return self.df["date_ord"].to_numpy()


def save_metadata(path: Path, payload: dict) -> None:
    with path.open("w") as f:
        json.dump(payload, f, indent=2)


def load_trajectory(path: str | Path, name: str) -> Trajectory:
    df = pd.read_parquet(path)
    return Trajectory(name=name, df=df)


def build_patents_trajectory(
    patents_aug_dir: str | Path,
    output_path: str | Path,
    start_year: int = 1980,
    end_year: int = 2020,
) -> Path:
    output_path = Path(output_path)
    ensure_dir(output_path.parent)

    print(f"[patents] building trajectory -> {output_path}")
    print(f"[patents] source directory: {patents_aug_dir}")
    frames = []
    n_files = 0
    for path in sorted(Path(patents_aug_dir).glob("*.csv")):
        if not path.stem.isdigit():
            continue
        file_year = int(path.stem)
        if file_year < start_year or file_year > end_year:
            continue
        df = pd.read_csv(
            path,
            usecols=[
                "publication.date",
                "number of patents",
                "unique ipc combinations",
                "unique authors (only first)",
            ],
        )
        df["publication.date"] = pd.to_datetime(df["publication.date"], errors="coerce")
        df = df.dropna(subset=["publication.date"]).copy()
        frames.append(df)
        n_files += 1
        print(f"[patents] loaded {path.name}: {len(df):,} rows")

    print(f"[patents] concatenating {n_files} yearly files")
    full = pd.concat(frames, ignore_index=True)
    print(f"[patents] total rows before global sort: {len(full):,}")
    full = full.sort_values("publication.date", kind="stable").reset_index(drop=True)
    print(f"[patents] global sort completed")
    full["date_ord"] = full["publication.date"].map(pd.Timestamp.toordinal).astype(np.int32)
    full["tau"] = full["publication.date"].map(decimal_year).astype(np.float32)
    full["year"] = full["publication.date"].dt.year.astype(np.uint16)
    full["t"] = np.arange(1, len(full) + 1, dtype=np.int64)
    full = full.rename(
        columns={
            "publication.date": "date",
            "number of patents": "t_raw",
            "unique ipc combinations": "D",
            "unique authors (only first)": "Dw",
        }
    )
    full["D"] = full["D"].astype(np.uint32)
    full["Dw"] = full["Dw"].astype(np.uint32)
    full = full[["date", "date_ord", "tau", "year", "t", "D", "Dw"]]
    print(f"[patents] writing parquet: {len(full):,} events")
    full.to_parquet(output_path, index=False)
    print(f"[patents] done")

    save_metadata(
        metadata_path(output_path.parent, "patents"),
        {
            "source": str(patents_aug_dir),
            "n_events": int(len(full)),
            "start_year": start_year,
            "end_year": end_year,
        },
    )
    return output_path


def _paper_is_valid(
    date_value: object,
    authors: object,
    keywords: object,
    topics: object,
    require_keywords: bool = True,
) -> bool:
    if pd.isna(date_value):
        return False
    if not isinstance(authors, (list, tuple, np.ndarray)) or len(authors) == 0:
        return False
    if not isinstance(topics, (list, tuple, np.ndarray)) or len(topics) == 0:
        return False
    first_author = str(authors[0]).strip()
    if not first_author:
        return False
    if require_keywords:
        if not isinstance(keywords, (list, tuple, np.ndarray)) or len(keywords) == 0:
            return False
        combo = normalize_keyword_combo(keywords)
        if not combo:
            return False
    return True


def _paper_year_mask(date_series: pd.Series, start_year: int, end_year: int) -> pd.Series:
    years = pd.to_numeric(date_series.astype("string").str.slice(0, 4), errors="coerce")
    return years.between(start_year, end_year)


def _chunked_lookup_map(
    conn: sqlite3.Connection,
    table_name: str,
    key_column: str,
    values: list[str],
    chunk_size: int = 500,
) -> dict[str, int]:
    result: dict[str, int] = {}
    for start in range(0, len(values), chunk_size):
        chunk = values[start : start + chunk_size]
        placeholders = ",".join("?" for _ in chunk)
        query = (
            f"SELECT {key_column}, id FROM {table_name} "
            f"WHERE {key_column} IN ({placeholders})"
        )
        for key, row_id in conn.execute(query, chunk):
            result[key] = row_id
    return result


def _init_papers_db(conn: sqlite3.Connection) -> None:
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.execute("PRAGMA temp_store=MEMORY;")
    conn.execute("PRAGMA cache_size=-200000;")
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS authors (
            id INTEGER PRIMARY KEY,
            author TEXT UNIQUE
        );
        CREATE TABLE IF NOT EXISTS combos (
            id INTEGER PRIMARY KEY,
            combo TEXT UNIQUE
        );
        CREATE TABLE IF NOT EXISTS events (
            event_id INTEGER PRIMARY KEY AUTOINCREMENT,
            date_text TEXT,
            date_ord INTEGER,
            tau REAL,
            author_id INTEGER,
            combo_id INTEGER
        );
        """
    )
    conn.commit()


def build_papers_trajectory(
    papers_dir: str | Path,
    output_path: str | Path,
    staging_dir: str | Path,
    start_year: int = 1980,
    end_year: int = 2020,
    subfields: list[int] | None = None,
    cleanup_stage: bool = True,
    require_keywords: bool = True,
) -> Path:
    output_path = Path(output_path)
    staging_dir = ensure_dir(staging_dir)
    ensure_dir(output_path.parent)

    mode_suffix = papers_mode_suffix(require_keywords)
    db_path = staging_dir / f"papers_trajectory_build_{mode_suffix}.sqlite"

    # Always rebuild from source. Both the parquet and the staging database of
    # any previous run are removed first: a leftover database would otherwise
    # be picked up as already-staged events.
    if output_path.exists():
        output_path.unlink()
    if db_path.exists():
        db_path.unlink()

    print(f"[papers] building trajectory -> {output_path}")
    print(f"[papers] mode: {mode_suffix}")
    print(f"[papers] source directory: {papers_dir}")
    print(f"[papers] staging database: {db_path}")
    conn = sqlite3.connect(db_path)
    _init_papers_db(conn)

    papers_dir = Path(papers_dir)
    if subfields is None:
        subfields = OPENALEX_SUBFIELDS
    parquet_files = [papers_dir / f"merged_df_subfield_{subfield}.parquet" for subfield in subfields]
    parquet_files = [path for path in parquet_files if path.exists()]
    print(f"[papers] parquet files to scan: {len(parquet_files)}")
    print(f"[papers] selected subfields: {sorted(subfields)}")

    total_inserted = 0
    total_valid = 0
    total_scanned = 0
    for file_idx, path in enumerate(parquet_files, start=1):
        print(f"[papers] [{file_idx}/{len(parquet_files)}] reading {path.name}")
        df = pd.read_parquet(path, columns=["date", "authors", "keywords", "topics"])
        total_scanned += len(df)
        mask_year = _paper_year_mask(df["date"], start_year, end_year)
        if not mask_year.any():
            print(f"[papers] [{file_idx}/{len(parquet_files)}] no rows in {start_year}-{end_year}, skip")
            continue
        df = df.loc[mask_year, ["date", "authors", "keywords", "topics"]].reset_index(drop=True)
        print(f"[papers] [{file_idx}/{len(parquet_files)}] rows in window: {len(df):,}")

        valid_rows = []
        unique_authors = set()
        unique_combos = set()
        for date_value, authors, keywords, topics in zip(
            df["date"], df["authors"], df["keywords"], df["topics"]
        ):
            if not _paper_is_valid(
                date_value,
                authors,
                keywords,
                topics,
                require_keywords=require_keywords,
            ):
                continue
            ts = pd.to_datetime(date_value, errors="coerce")
            if pd.isna(ts):
                continue
            first_author = str(authors[0]).strip()
            combo = None
            if isinstance(keywords, (list, tuple, np.ndarray)) and len(keywords) > 0:
                combo = normalize_keyword_combo(keywords)
                if not combo:
                    combo = None
            date_text = ts.strftime("%Y-%m-%d")
            valid_rows.append(
                (
                    date_text,
                    int(ts.toordinal()),
                    float(decimal_year(ts)),
                    first_author,
                    combo,
                )
            )
            unique_authors.add(first_author)
            if combo is not None:
                unique_combos.add(combo)

        if not valid_rows:
            print(f"[papers] [{file_idx}/{len(parquet_files)}] no valid rows after filtering")
            continue

        total_valid += len(valid_rows)
        print(
            f"[papers] [{file_idx}/{len(parquet_files)}] valid rows: {len(valid_rows):,} | "
            f"unique authors chunk: {len(unique_authors):,} | unique combos chunk: {len(unique_combos):,}"
        )

        conn.executemany(
            "INSERT OR IGNORE INTO authors(author) VALUES (?)",
            [(author,) for author in unique_authors],
        )
        conn.executemany(
            "INSERT OR IGNORE INTO combos(combo) VALUES (?)",
            [(combo,) for combo in unique_combos],
        )
        conn.commit()

        author_map = _chunked_lookup_map(conn, "authors", "author", sorted(unique_authors))
        combo_map = _chunked_lookup_map(conn, "combos", "combo", sorted(unique_combos))

        event_rows = [
            (
                date_text,
                date_ord,
                tau,
                author_map[first_author],
                combo_map[combo] if combo is not None else None,
            )
            for date_text, date_ord, tau, first_author, combo in valid_rows
        ]
        conn.executemany(
            "INSERT INTO events(date_text, date_ord, tau, author_id, combo_id) VALUES (?, ?, ?, ?, ?)",
            event_rows,
        )
        conn.commit()
        total_inserted += len(event_rows)
        n_authors_db = conn.execute("SELECT COUNT(*) FROM authors").fetchone()[0]
        n_combos_db = conn.execute("SELECT COUNT(*) FROM combos").fetchone()[0]
        print(
            f"[papers] [{file_idx}/{len(parquet_files)}] staged events so far: {total_inserted:,} | "
            f"authors db: {n_authors_db:,} | combos db: {n_combos_db:,}"
        )

    print(f"[papers] raw rows scanned: {total_scanned:,}")
    print(f"[papers] valid rows staged: {total_valid:,}")
    print(f"[papers] creating index for ordered scan")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_events_date ON events(date_ord, event_id);")
    conn.commit()
    n_events = conn.execute("SELECT COUNT(*) FROM events").fetchone()[0]
    print(f"[papers] events in staging db: {n_events:,}")

    query = (
        "SELECT date_text, date_ord, tau, author_id, combo_id "
        "FROM events ORDER BY date_ord, event_id"
    )
    cursor = conn.execute(query)
    print(f"[papers] ordered cursor ready, writing final parquet")

    schema = pa.schema(
        [
            ("date", pa.timestamp("ns")),
            ("date_ord", pa.int32()),
            ("tau", pa.float32()),
            ("year", pa.uint16()),
            ("t", pa.int64()),
            ("D", pa.uint32()),
            ("Dw", pa.uint32()),
        ]
    )
    writer = pq.ParquetWriter(output_path, schema=schema, compression="zstd")

    seen_authors: set[int] = set()
    seen_combos: set[int] = set()
    d_count = 0
    w_count = 0
    t_count = 0
    batch_size = 200_000
    progress_every = 1_000_000

    batch_date = []
    batch_date_ord = []
    batch_tau = []
    batch_year = []
    batch_t = []
    batch_D = []
    batch_Dw = []

    for date_text, date_ord, tau, author_id, combo_id in cursor:
        t_count += 1
        if combo_id is not None and combo_id not in seen_combos:
            seen_combos.add(combo_id)
            d_count += 1
        if author_id not in seen_authors:
            seen_authors.add(author_id)
            w_count += 1

        batch_date.append(np.datetime64(date_text, "ns"))
        batch_date_ord.append(int(date_ord))
        batch_tau.append(float(tau))
        batch_year.append(int(date_text[:4]))
        batch_t.append(t_count)
        batch_D.append(d_count)
        batch_Dw.append(w_count)

        if len(batch_t) >= batch_size:
            table = pa.table(
                {
                    "date": pa.array(batch_date, type=pa.timestamp("ns")),
                    "date_ord": pa.array(batch_date_ord, type=pa.int32()),
                    "tau": pa.array(batch_tau, type=pa.float32()),
                    "year": pa.array(batch_year, type=pa.uint16()),
                    "t": pa.array(batch_t, type=pa.int64()),
                    "D": pa.array(batch_D, type=pa.uint32()),
                    "Dw": pa.array(batch_Dw, type=pa.uint32()),
                },
                schema=schema,
            )
            writer.write_table(table)
            batch_date.clear()
            batch_date_ord.clear()
            batch_tau.clear()
            batch_year.clear()
            batch_t.clear()
            batch_D.clear()
            batch_Dw.clear()
        if t_count % progress_every == 0:
            print(
                f"[papers] final pass progress: {t_count:,}/{n_events:,} events | "
                f"D={d_count:,} | Dw={w_count:,}"
            )

    if batch_t:
        table = pa.table(
            {
                "date": pa.array(batch_date, type=pa.timestamp("ns")),
                "date_ord": pa.array(batch_date_ord, type=pa.int32()),
                "tau": pa.array(batch_tau, type=pa.float32()),
                "year": pa.array(batch_year, type=pa.uint16()),
                "t": pa.array(batch_t, type=pa.int64()),
                "D": pa.array(batch_D, type=pa.uint32()),
                "Dw": pa.array(batch_Dw, type=pa.uint32()),
            },
            schema=schema,
        )
        writer.write_table(table)
    print(f"[papers] final pass completed: {t_count:,} events")

    writer.close()
    conn.close()
    print(f"[papers] parquet written: {output_path}")

    save_metadata(
        papers_metadata_path(output_path.parent, require_keywords),
        {
            "source": str(papers_dir),
            "n_events": int(t_count),
            "start_year": start_year,
            "end_year": end_year,
            "staging_db": str(db_path),
            "mode": mode_suffix,
            "require_keywords": require_keywords,
        },
    )

    if cleanup_stage and db_path.exists():
        print(f"[papers] removing staging db")
        db_path.unlink()

    print(f"[papers] done")
    return output_path


def build_all_trajectories(
    trajectory_dir: str | Path,
    patents_aug_dir: str | Path,
    papers_dir: str | Path,
    staging_dir: str | Path,
    papers_subfields: list[int] | None = None,
    cleanup_stage: bool = True,
    filter_empty_keywords: bool = True,
) -> dict[str, Path]:
    paths = trajectory_paths(trajectory_dir, filter_empty_keywords=filter_empty_keywords)
    build_patents_trajectory(patents_aug_dir, paths["patents"])
    build_papers_trajectory(
        papers_dir,
        paths["papers"],
        staging_dir=staging_dir,
        subfields=papers_subfields,
        cleanup_stage=cleanup_stage,
        require_keywords=filter_empty_keywords,
    )
    return paths


def load_all_trajectories(
    trajectory_dir: str | Path,
    filter_empty_keywords: bool = True,
) -> dict[str, Trajectory]:
    paths = trajectory_paths(trajectory_dir, filter_empty_keywords=filter_empty_keywords)
    return {
        "patents": load_trajectory(paths["patents"], "patents"),
        "papers": load_trajectory(paths["papers"], f"papers[filter_empty_keywords={filter_empty_keywords}]"),
    }


def compress_by_date(traj: Trajectory) -> pd.DataFrame:
    df = traj.df
    mask = np.r_[df["date_ord"].to_numpy()[1:] != df["date_ord"].to_numpy()[:-1], True]
    return df.loc[mask, ["t", "tau", "D", "Dw", "year"]].reset_index(drop=True)


def annual_summary(traj: Trajectory) -> pd.DataFrame:
    df = traj.df
    grouped = df.groupby("year", sort=True)
    out = grouped.agg(
        tau=("year", "first"),
        delta_t=("t", "size"),
        D=("D", "last"),
        Dw=("Dw", "last"),
        t=("t", "last"),
    )
    return out.reset_index()


def display_stride(n: int, max_points: int = 50_000) -> int:
    return max(1, int(math.ceil(n / max_points)))


def beta_eff_binned(
    t: np.ndarray,
    D: np.ndarray,
    n_log_bins: int = 400,
    smooth_window: int = 11,
    trim_windows: int = 1,
    t_min: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Local effective exponent on the log-time mesh used by :func:`beta_eff`.

    Returns the estimate on its own mesh instead of interpolating it back onto
    the event-level trajectory, so that the boundary bins can be discarded.
    ``trim_windows`` drops that many smoothing windows at each end: the local
    slope is biased there because the smoothing pads with edge values and the
    gradient uses one-sided differences. ``t_min`` additionally discards the
    initial transient, where ``D(t)`` is still dominated by the first events.
    """

    t = np.asarray(t, dtype=np.float64)
    D = np.asarray(D, dtype=np.float64)
    mask = (t > 1) & (D > 0)
    if not np.any(mask):
        return np.empty(0), np.empty(0)

    t_valid = t[mask]
    D_valid = D[mask]
    change_mask = np.r_[True, D_valid[1:] != D_valid[:-1]]
    t_change = t_valid[change_mask]
    D_change = D_valid[change_mask]
    if len(t_change) < 3:
        return np.empty(0), np.empty(0)

    log_t_change = np.log(t_change)
    log_edges = np.linspace(log_t_change.min(), log_t_change.max(), n_log_bins + 1)
    bin_ids = np.digitize(log_t_change, log_edges[1:-1], right=False)

    t_binned = []
    D_binned = []
    for bin_id in range(n_log_bins):
        idx = np.where(bin_ids == bin_id)[0]
        if idx.size == 0:
            continue
        t_binned.append(t_change[idx[-1]])
        D_binned.append(D_change[idx[-1]])

    t_binned = np.asarray(t_binned, dtype=np.float64)
    D_binned = np.asarray(D_binned, dtype=np.float64)
    if len(t_binned) < 3:
        return np.empty(0), np.empty(0)

    log_t = np.log(t_binned)
    log_D = np.log(D_binned)
    if smooth_window >= 3 and len(log_D) >= smooth_window:
        smooth_window = smooth_window if smooth_window % 2 == 1 else smooth_window + 1
        kernel = np.ones(smooth_window, dtype=np.float64) / smooth_window
        pad = smooth_window // 2
        log_D = np.convolve(np.pad(log_D, pad, mode="edge"), kernel, mode="valid")

    beta_binned = np.gradient(log_D, log_t)

    trim = max(int(trim_windows) * int(smooth_window), 0)
    if trim > 0 and len(t_binned) > 2 * trim + 2:
        t_binned = t_binned[trim:-trim]
        beta_binned = beta_binned[trim:-trim]
    if t_min is not None:
        keep = t_binned >= float(t_min)
        t_binned = t_binned[keep]
        beta_binned = beta_binned[keep]
    return t_binned, beta_binned


def regression_metrics(y: np.ndarray, yhat: np.ndarray, k: int) -> dict[str, float]:
    y = np.asarray(y, dtype=np.float64)
    yhat = np.asarray(yhat, dtype=np.float64)
    residuals = y - yhat
    sse = float(np.sum(residuals**2))
    n = int(len(y))
    mse = sse / n
    rmse = math.sqrt(mse)
    sst = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - sse / sst if sst > 0 else np.nan
    sigma2 = max(sse / n, np.finfo(float).eps)
    aic = n * math.log(sigma2) + 2 * k
    bic = n * math.log(sigma2) + k * math.log(n)
    return {"sse": sse, "mse": mse, "rmse": rmse, "r2": r2, "aic": aic, "bic": bic}


def log_regression_metrics(y: np.ndarray, yhat: np.ndarray, k: int) -> dict[str, float]:
    y = np.asarray(y, dtype=np.float64)
    yhat = np.asarray(yhat, dtype=np.float64)
    mask = (y > 0) & (yhat > 0)
    return regression_metrics(np.log(y[mask]), np.log(yhat[mask]), k)


def _through_origin_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    beta = float(np.dot(x, y) / np.dot(x, x))
    resid = y - beta * x
    n = len(x)
    sigma2 = float(np.sum(resid**2) / max(n - 1, 1))
    se = math.sqrt(sigma2 / np.dot(x, x))
    t_stat = beta / se if se > 0 else np.inf
    pvalue = 2 * stats.t.sf(abs(t_stat), df=max(n - 1, 1))
    return beta, se, float(pvalue)


def fit_intrinsic_models(t: np.ndarray, D: np.ndarray) -> pd.DataFrame:
    t = np.asarray(t, dtype=np.float64)
    D = np.asarray(D, dtype=np.float64)
    mask = t > 1
    t = t[mask]
    D = D[mask]

    rows = []

    x_linear = t
    a_linear, se_a_linear, p_linear = _through_origin_fit(x_linear, D)
    yhat_linear = a_linear * x_linear
    rows.append(
        {
            "model": "a * t",
            "a": a_linear,
            "se_a": se_a_linear,
            "beta": np.nan,
            "se_beta": np.nan,
            "p_value_main": p_linear,
            **regression_metrics(D, yhat_linear, k=1),
            **{f"log_{k}": v for k, v in log_regression_metrics(D, yhat_linear, k=1).items()},
        }
    )

    x_marginal = t / np.log(t)
    a_marginal, se_a_marginal, p_marginal = _through_origin_fit(x_marginal, D)
    yhat_marginal = a_marginal * x_marginal
    rows.append(
        {
            "model": "a * t / log(t)",
            "a": a_marginal,
            "se_a": se_a_marginal,
            "beta": np.nan,
            "se_beta": np.nan,
            "p_value_main": p_marginal,
            **regression_metrics(D, yhat_marginal, k=1),
            **{f"log_{k}": v for k, v in log_regression_metrics(D, yhat_marginal, k=1).items()},
        }
    )

    log_t = np.log(t)
    log_D = np.log(D)
    ols = sm.OLS(log_D, sm.add_constant(log_t)).fit()
    beta_power = float(ols.params[1])
    a_power = float(math.exp(ols.params[0]))
    # Delta method: a = exp(c) implies se_a = a * se_c.
    se_a_power = a_power * float(ols.bse[0])
    yhat_power = a_power * np.power(t, beta_power)
    rows.append(
        {
            "model": "a * t^beta",
            "a": a_power,
            "se_a": se_a_power,
            "beta": beta_power,
            "se_beta": float(ols.bse[1]),
            "p_value_main": float(ols.pvalues[1]),
            **regression_metrics(D, yhat_power, k=2),
            **{f"log_{k}": v for k, v in log_regression_metrics(D, yhat_power, k=2).items()},
        }
    )

    return pd.DataFrame(rows).sort_values("aic").reset_index(drop=True)


def predict_intrinsic_model(model_name: str, params: pd.Series, t: np.ndarray) -> np.ndarray:
    t = np.asarray(t, dtype=np.float64)
    if model_name == "a * t":
        return params["a"] * t
    if model_name == "a * t / log(t)":
        return params["a"] * t / np.log(t)
    if model_name == "a * t^beta":
        return params["a"] * np.power(t, params["beta"])
    raise ValueError(f"Unknown model: {model_name}")


def fit_exponential_vs_tau(
    tau: np.ndarray,
    series: np.ndarray,
    tau_min: float | None = None,
    tau_max: float | None = None,
) -> dict[str, float]:
    tau = np.asarray(tau, dtype=np.float64)
    series = np.asarray(series, dtype=np.float64)
    mask = series > 0
    if tau_min is not None:
        mask &= tau >= tau_min
    if tau_max is not None:
        mask &= tau <= tau_max
    X = sm.add_constant(tau[mask])
    fit = sm.OLS(np.log(series[mask]), X).fit()
    slope = float(fit.params[1])
    intercept = float(fit.params[0])
    yhat = np.exp(intercept + slope * tau[mask])
    metrics = regression_metrics(series[mask], yhat, k=2)
    return {
        "intercept": intercept,
        "slope": slope,
        "p_value_slope": float(fit.pvalues[1]),
        "tau_min": float(np.min(tau[mask])),
        "tau_max": float(np.max(tau[mask])),
        **metrics,
    }


def branching_models(D: np.ndarray, Dw: np.ndarray) -> pd.DataFrame:
    x = np.asarray(D, dtype=np.float64)
    y = np.asarray(Dw, dtype=np.float64)
    mask = (x > 1) & (y > 0)
    x = x[mask]
    y = y[mask]

    rows = []

    b_lin, se_b_lin, p_lin = _through_origin_fit(x, y)
    yhat_lin = b_lin * x
    rows.append(
        {
            "model": "linear: b * D",
            "a": np.nan,
            "se_a": np.nan,
            "b": b_lin,
            "se_b": se_b_lin,
            "beta": np.nan,
            "se_beta": np.nan,
            "p_value_main": p_lin,
            **regression_metrics(y, yhat_lin, k=1),
            **{f"log_{k}": v for k, v in log_regression_metrics(y, yhat_lin, k=1).items()},
        }
    )

    fit_pow = sm.OLS(np.log(y), sm.add_constant(np.log(x))).fit()
    beta_pow = float(fit_pow.params[1])
    a_pow = float(math.exp(fit_pow.params[0]))
    yhat_pow = a_pow * np.power(x, beta_pow)
    rows.append(
        {
            "model": "power: a * D^beta",
            "a": a_pow,
            "se_a": a_pow * float(fit_pow.bse[0]),
            "b": np.nan,
            "se_b": np.nan,
            "beta": beta_pow,
            "se_beta": float(fit_pow.bse[1]),
            "p_value_main": float(fit_pow.pvalues[1]),
            **regression_metrics(y, yhat_pow, k=2),
            **{f"log_{k}": v for k, v in log_regression_metrics(y, yhat_pow, k=2).items()},
        }
    )

    fit_log = sm.OLS(y, sm.add_constant(np.log(x))).fit()
    a_log = float(fit_log.params[0])
    b_log = float(fit_log.params[1])
    yhat_log = a_log + b_log * np.log(x)
    rows.append(
        {
            "model": "log: a + b * log(D)",
            "a": a_log,
            "se_a": float(fit_log.bse[0]),
            "b": b_log,
            "se_b": float(fit_log.bse[1]),
            "beta": np.nan,
            "se_beta": np.nan,
            "p_value_main": float(fit_log.pvalues[1]),
            **regression_metrics(y, yhat_log, k=2),
            **{f"log_{k}": v for k, v in log_regression_metrics(y, yhat_log, k=2).items()},
        }
    )

    x_scaled = x / np.max(x)
    fit_exp = sm.OLS(np.log(y), sm.add_constant(x_scaled)).fit()
    a_exp = float(math.exp(fit_exp.params[0]))
    b_exp = float(fit_exp.params[1])
    yhat_exp = a_exp * np.exp(b_exp * x_scaled)
    rows.append(
        {
            "model": "exp: a * exp(b * D_scaled)",
            "a": a_exp,
            "se_a": a_exp * float(fit_exp.bse[0]),
            "b": b_exp,
            "se_b": float(fit_exp.bse[1]),
            "beta": np.nan,
            "se_beta": np.nan,
            "p_value_main": float(fit_exp.pvalues[1]),
            **regression_metrics(y, yhat_exp, k=2),
            **{f"log_{k}": v for k, v in log_regression_metrics(y, yhat_exp, k=2).items()},
        }
    )

    return pd.DataFrame(rows).sort_values("aic").reset_index(drop=True)


def predict_branching_model(model_name: str, params: pd.Series, D: np.ndarray) -> np.ndarray:
    D = np.asarray(D, dtype=np.float64)
    if model_name == "linear: b * D":
        return params["b"] * D
    if model_name == "power: a * D^beta":
        return params["a"] * np.power(D, params["beta"])
    if model_name == "log: a + b * log(D)":
        return params["a"] + params["b"] * np.log(D)
    if model_name == "exp: a * exp(b * D_scaled)":
        scaled = D / np.max(D)
        return params["a"] * np.exp(params["b"] * scaled)
    raise ValueError(f"Unknown model: {model_name}")


def block_b_summary(
    D: np.ndarray,
    Dw: np.ndarray,
    tau: np.ndarray,
    block_size: int,
    smooth_blocks: int = 9,
) -> pd.DataFrame:
    D = np.asarray(D, dtype=np.float64)
    Dw = np.asarray(Dw, dtype=np.float64)
    tau = np.asarray(tau, dtype=np.float64)

    rows = []
    n = len(D)
    for start in range(0, n - block_size + 1, block_size):
        stop = start + block_size
        dD = D[stop - 1] - D[start]
        dDw = Dw[stop - 1] - Dw[start]
        if dD <= 0:
            continue
        rows.append(
            {
                "tau_center": float(np.mean(tau[start:stop])),
                "b_local": float(dDw / dD),
                "delta_D": float(dD),
                "delta_Dw": float(dDw),
                "start_t": int(start + 1),
                "end_t": int(stop),
            }
        )

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    out["Dw_over_D"] = out["delta_Dw"] / out["delta_D"]

    window = max(1, int(smooth_blocks))
    if window % 2 == 0:
        window += 1

    out["b_mean"] = out["b_local"].rolling(window, center=True, min_periods=1).mean()
    out["b_std"] = out["b_local"].rolling(window, center=True, min_periods=1).std(ddof=0)
    out["b_q16"] = out["b_local"].rolling(window, center=True, min_periods=1).quantile(0.16)
    out["b_q84"] = out["b_local"].rolling(window, center=True, min_periods=1).quantile(0.84)

    out["ratio_mean"] = out["Dw_over_D"].rolling(window, center=True, min_periods=1).mean()
    out["ratio_std"] = out["Dw_over_D"].rolling(window, center=True, min_periods=1).std(ddof=0)
    out["ratio_q16"] = out["Dw_over_D"].rolling(window, center=True, min_periods=1).quantile(0.16)
    out["ratio_q84"] = out["Dw_over_D"].rolling(window, center=True, min_periods=1).quantile(0.84)

    return out
