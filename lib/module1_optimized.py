from __future__ import annotations

import json
import pickle
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from analysis_utils import _load_papers_global_frame, normalize_tokens
from utils import OPENALEX_SUBFIELDS, ensure_dir, normalize_keyword_combo, reduction_indices


def _compute_intervals(author_dates: dict[str, list[dict]]) -> dict[str, list[dict]]:
    def compute_intervals_with_dates(events_list: list[dict], author_name: str) -> list[dict]:
        if len(events_list) < 2:
            return []
        dates = sorted(pd.to_datetime(ev["date"]) for ev in events_list)
        intervals = []
        for d1, d2 in zip(dates[:-1], dates[1:]):
            dt = (d2 - d1).days / 365.25
            if 0 < dt < 50:
                intervals.append(
                    {
                        "interval_years": dt,
                        "start_date": d1,
                        "end_date": d2,
                        "author": author_name,
                    }
                )
        return intervals

    all_intervals = []
    global_only_intervals = []
    novelties_intervals = []
    no_novelty_intervals = []
    for author, events in author_dates.items():
        all_intervals.extend(compute_intervals_with_dates(events, author))
        global_only_intervals.extend(
            compute_intervals_with_dates(
                [ev for ev in events if ev["novelty_type"] == "global_novelty"],
                author,
            )
        )
        novelties_intervals.extend(
            compute_intervals_with_dates(
                [ev for ev in events if ev["novelty_type"] in ("global_novelty", "individual_novelty")],
                author,
            )
        )
        no_novelty_intervals.extend(
            compute_intervals_with_dates(
                [ev for ev in events if ev["novelty_type"] == "no_novelty"],
                author,
            )
        )
    return {
        "all_intervals": all_intervals,
        "global_only_intervals": global_only_intervals,
        "novelties_intervals": novelties_intervals,
        "no_novelty_intervals": no_novelty_intervals,
    }


def _write_pickle(path: Path, payload: object) -> None:
    with path.open("wb") as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)


def build_papers_module1_outputs(
    papers_dir: str | Path,
    output_dir: str | Path,
    subfields: list[int] | None = None,
    start_year: int = 1980,
    end_year: int = 2020,
    filter_empty_keywords: bool = True,
    max_papers_per_year: int = 100,
    n_points: int = 10_000,
) -> dict[str, object]:
    print("[module1/papers] building outputs")
    papers_dir = Path(papers_dir)
    output_dir = ensure_dir(output_dir)
    if subfields is None:
        subfields = OPENALEX_SUBFIELDS
    print("[module1/papers] loading globally ordered papers frame")
    df = _load_papers_global_frame(
        papers_dir=papers_dir,
        subfields=subfields,
        start_year=start_year,
        end_year=end_year,
    )
    print(f"[module1/papers] global frame loaded: {len(df):,} rows before filters")

    valid_topics = df["topics"].apply(lambda x: isinstance(x, (list, tuple, np.ndarray)) and len(x) > 0)
    valid_authors = df["authors"].apply(lambda x: isinstance(x, (list, tuple, np.ndarray)) and len(x) > 0)
    valid_keywords = df["keywords"].apply(lambda x: isinstance(x, (list, tuple, np.ndarray)) and len(x) > 0)
    keep_mask = valid_topics & valid_authors
    if filter_empty_keywords:
        keep_mask &= valid_keywords
    df = df.loc[keep_mask].copy()
    df["year"] = df["date"].dt.year.astype(int)

    print("[module1/papers] computing excluded prolific first authors")
    df["first_author"] = df["authors"].apply(
        lambda a: str(a[0]).strip() if isinstance(a, (list, tuple, np.ndarray)) and len(a) > 0 else ""
    )
    counts = (
        df.loc[df["first_author"] != "", ["year", "first_author"]]
        .groupby(["year", "first_author"])
        .size()
    )
    excluded_authors = set(counts[counts > max_papers_per_year].reset_index()["first_author"])
    excluded_count = len(excluded_authors)
    total_before = len(df)
    print(
        f"[module1/papers] total rows before exclusion: {total_before:,} | "
        f"excluded authors: {excluded_count:,}"
    )
    if excluded_authors:
        df = df.loc[~df["first_author"].isin(excluded_authors)].copy()
    total_after = len(df)
    print(f"[module1/papers] total rows after exclusion: {total_after:,}")

    log_indices, lin_indices = reduction_indices(int(total_after), n_points=n_points)
    next_log = 0
    next_lin = 0
    seen_combo: set[str] = set()
    seen_author: set[str] = set()
    seen_keywords: set[str] = set()
    seen_topic_global: set[str] = set()
    seen_topic_individual: dict[str, set[str]] = defaultdict(set)
    combo_counter: Counter[str] = Counter()
    author_counter: Counter[str] = Counter()
    author_dates: dict[str, list[dict]] = defaultdict(list)
    years_log = []
    years_lin = []
    combos_log = []
    combos_lin = []
    authors_log = []
    authors_lin = []
    rows_log = []
    rows_lin = []
    delta_rows: Counter[int] = Counter()
    authors_last_per_year: dict[int, int] = {}

    event_idx = -1
    progress_step = 500_000
    for row in df.itertuples(index=False):
        event_idx += 1
        if event_idx == 0 or (event_idx + 1) % progress_step == 0:
            print(
                f"[module1/papers] final pass {event_idx + 1:,}/{total_after:,} | "
                f"seen combos={len(seen_combo):,} | seen authors={len(seen_author):,}"
            )
        date_value = pd.Timestamp(row.date)
        year = int(row.year)
        first_author = row.first_author
        keywords = normalize_tokens(row.keywords)
        topics = normalize_tokens(row.topics)
        keyword_combo = normalize_keyword_combo(keywords) if keywords else ""
        topics_combo = normalize_keyword_combo(topics) if topics else ""

        delta_rows[year] += 1
        combo_counter[keyword_combo] += 1
        author_counter[first_author] += 1
        if keyword_combo:
            seen_combo.add(keyword_combo)
            for token in keyword_combo.split("\x1f"):
                if token:
                    seen_keywords.add(token)
        seen_author.add(first_author)
        authors_last_per_year[year] = len(seen_author)

        if next_log < len(log_indices) and event_idx == int(log_indices[next_log]):
            years_log.append(year)
            combos_log.append(len(seen_combo))
            authors_log.append(len(seen_author))
            rows_log.append(event_idx + 1)
            next_log += 1
        if next_lin < len(lin_indices) and event_idx == int(lin_indices[next_lin]):
            years_lin.append(year)
            combos_lin.append(len(seen_combo))
            authors_lin.append(len(seen_author))
            rows_lin.append(event_idx + 1)
            next_lin += 1

        if topics_combo not in seen_topic_global:
            novelty_type = "global_novelty"
            seen_topic_global.add(topics_combo)
            seen_topic_individual[first_author].add(topics_combo)
        elif topics_combo not in seen_topic_individual[first_author]:
            novelty_type = "individual_novelty"
            seen_topic_individual[first_author].add(topics_combo)
        else:
            novelty_type = "no_novelty"
        author_dates[first_author].append(
            {
                "date": date_value,
                "year": year,
                "topics": topics_combo.split("\x1f") if topics_combo else [],
                "novelty_type": novelty_type,
            }
        )
    print("[module1/papers] writing csv/pkl outputs")

    years_sorted = sorted(delta_rows)
    delta_rows_df = pd.DataFrame(
        {
            "delta_rows": [int(delta_rows[y]) for y in years_sorted],
            "num_first_authors": [int(authors_last_per_year.get(y, 0)) for y in years_sorted],
        }
    )
    delta_rows_df.to_csv(output_dir / "delta_rows_papers.csv", index=False)

    df_combo_freq = pd.DataFrame(
        {
            "rank": np.arange(1, len(combo_counter) + 1),
            "combo": [combo.split("\x1f") if combo else [] for combo, _ in combo_counter.most_common()],
            "frequency": [count for _, count in combo_counter.most_common()],
        }
    )
    df_combo_freq.to_csv(output_dir / "openalex_combo_freq.csv", index=False)

    df_author_freq = pd.DataFrame(
        {
            "rank": np.arange(1, len(author_counter) + 1),
            "author": [author for author, _ in author_counter.most_common()],
            "frequency": [count for _, count in author_counter.most_common()],
        }
    )
    df_author_freq.to_csv(output_dir / "openalex_author_freq.csv", index=False)

    pd.DataFrame(
        {
            "years": years_log,
            "num_unique_combinations": combos_log,
            "num_unique_authors": authors_log,
            "row_indices": rows_log,
        }
    ).to_csv(output_dir / "reduced_openalex_log.csv", index=False)

    pd.DataFrame(
        {
            "years": years_lin,
            "num_unique_combinations": combos_lin,
            "num_unique_authors": authors_lin,
            "row_indices": rows_lin,
        }
    ).to_csv(output_dir / "reduced_openalex_lin.csv", index=False)

    author_dates_path = output_dir / f"openalex_author_dates_{start_year}_{end_year}.pkl"
    intervals_path = output_dir / "openalex_intervals.pkl"
    _write_pickle(author_dates_path, dict(author_dates))
    _write_pickle(intervals_path, _compute_intervals(dict(author_dates)))

    meta = {
        "strategy": "global_frame_single_pass",
        "staged_rows": None,
        "total_rows_before_exclusion": int(total_before),
        "total_rows_after_exclusion": int(total_after),
        "excluded_authors": int(excluded_count),
        "unique_keywords": int(len(seen_keywords)),
        "unique_keyword_combinations": int(len(seen_combo)),
        "unique_first_authors": int(len(seen_author)),
        "filter_empty_keywords": bool(filter_empty_keywords),
        "max_papers_per_year": int(max_papers_per_year),
    }
    (output_dir / "module_1_papers_meta.json").write_text(json.dumps(meta, indent=2))
    print("[module1/papers] done")
    return meta

