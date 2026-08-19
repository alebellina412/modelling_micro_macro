# Raw data

Put the two raw archives here. They are third-party data and are not
redistributed with this repository; everything in this directory except this
file is excluded from version control.

```
data/
├── patents_EPO/          one CSV per year, 1980.csv … 2020.csv
└── papers_openalex/      one parquet per subfield, merged_df_subfield_{ID}.parquet
```

**`patents_EPO/YYYY.csv`** — one row per patent, with at least:

| column | content |
|---|---|
| `publication.date` | `YYYYMMDD` |
| `inventor` | list-like string of dicts, each with a `name` |
| `ipc` | list-like string of IPC codes |

**`papers_openalex/merged_df_subfield_{ID}.parquet`** — one row per paper, with
at least:

| column | content |
|---|---|
| `date` | publication date |
| `authors` | list of author names, the first one is the focal contributor |
| `keywords` | list of keywords; their combination is the novelty |
| `topics` | list of topics |

The eighteen subfield IDs, and the order in which the files are read, are in
the parameter cell of `module_1_papers.ipynb`.

If the archives already live elsewhere on your machine, symlink them rather
than copying:

```bash
ln -s /path/to/patents_EPO      data/patents_EPO
ln -s /path/to/papers_openalex  data/papers_openalex
```
