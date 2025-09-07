# Mouse↔Human Ortholog Pairing & Scatter

Quick README for `combine_ortholog_scatter.py`. This script pairs **human** and **mouse** differential expression (DE) results via **1:1 orthologs**, exports a tidy table, and draws a cross‑species scatter plot (x = mouse log2FC, y = human log2FC) with a linear fit and Pearson **r**.

---

## What it does

* Reads **two CSVs** (human & mouse).
* Ensures each has: **first column = gene symbol**, plus **`log2FoldChange`** and **`padj`** columns (order doesn’t matter).
* Builds a **1:1 human→mouse** ortholog map via Ensembl BioMart *(or uses your own precomputed ortholog file).*
* Outputs a CSV with exactly these columns:

  * `gene_mouse, gene_human, log2FC_mouse, log2FC_human, padj_mouse, padj_human`
* Produces a PNG scatter filtered by significance (default **padj < 0.05 in both species**), with an OLS line and **Pearson r** annotation.

---

## Requirements

* Python 3.8+
* Packages: `pandas`, `requests`, `matplotlib`, `numpy`

```bash
pip install pandas requests matplotlib numpy
```

> **Note:** Online ortholog fetching requires outbound HTTPS to Ensembl BioMart. If your environment is offline, generate an ortholog CSV first (see below) and pass it via `--orthologs`.

---

## Inputs

* **Human CSV**: first column = human gene symbol; must contain `log2FoldChange` and `padj`.
* **Mouse CSV**: first column = mouse gene symbol; must contain `log2FoldChange` and `padj`.
* *(Optional)* **Ortholog CSV**: columns **`gene_human,gene_mouse`** (1:1 only). If omitted, the script queries Ensembl BioMart.

---

## Basic usage

```bash
python combine_ortholog_scatter.py \
  --human path/to/human.csv \
  --mouse path/to/mouse.csv \
  --out-csv output/paired_mouse_human.csv \
  --out-plot output/mouse_human_scatter.png
```

This will:

1. Fetch 1:1 human→mouse orthologs from BioMart.
2. Join human+mouse DE by orthologs.
3. Save `output/paired_mouse_human.csv` and `output/mouse_human_scatter.png`.

---

## Using a precomputed ortholog map (no network)

```bash
python combine_ortholog_scatter.py \
  --human path/to/human.csv \
  --mouse path/to/mouse.csv \
  --orthologs path/to/orthologs_1to1.csv \
  --out-csv output/paired_mouse_human.csv \
  --out-plot output/mouse_human_scatter.png
```

**Expected ortholog CSV columns:**

```
gene_human,gene_mouse
EPCAM,Epcam
KRT8,Krt8
...
```

---

## CLI options

```
--human <file>            Human DE CSV (required)
--mouse <file>            Mouse DE CSV (required)
--orthologs <file>        Optional 1:1 map (gene_human,gene_mouse). Skips BioMart.
--out-csv <file>          Output paired table (default: paired_mouse_human.csv)
--out-plot <file>         Output scatter PNG (default: mouse_human_scatter.png)
--padj-threshold <float>  Significance cutoff for plotting (default: 0.05)
--sig-mode <mode>         Which padj must pass for plotting: one of
                          {both|any|human|mouse} (default: both)
```

**Tip:** Use `--sig-mode any` to include genes significant in either species; `--sig-mode human` or `mouse` to gate by one side only.

---

## Output files

* **Paired table (CSV):**

  * Columns: `gene_mouse, gene_human, log2FC_mouse, log2FC_human, padj_mouse, padj_human`
  * Rows: genes present in both DE tables and mapped via 1:1 orthologs.
* **Scatter (PNG):**

  * X = mouse `log2FoldChange`, Y = human `log2FoldChange`.
  * Points filtered by `padj` per `--sig-mode` and `--padj-threshold`.
  * Linear regression line and **Pearson r** displayed.

---

## Quality & interpretation checks

* Ensure **contrast orientation** is consistent (e.g., epithelial/cyst vs comparator) so signs match biology.
* Prefer **shrunken log2FC** (e.g., DESeq2 `lfcShrink`) upstream to stabilize extremes.
* If the scatter looks sparse, try `--sig-mode any` or relax `--padj-threshold`.

---

## Troubleshooting

* **"No ortholog pairs found"**: Symbol mismatches or aliases. Consider passing a vetted ortholog CSV.
* **Empty plot (no points)**: Tight significance filter; relax with `--sig-mode any` or larger `--padj-threshold`.
* **BioMart errors/timeouts**: Use `--orthologs` to bypass network or retry later.

---

## Repro tips

* Commit the exact inputs and the produced `paired_mouse_human.csv` alongside the figure for full provenance.
* Record the CLI you used (copy from your shell history) in your analysis notes.

---

## Optional: building an ortholog CSV

If you don’t want the script to call BioMart, you can pre‑generate a strict 1:1 map (columns `gene_human,gene_mouse`) using your preferred pipeline or an external helper. Then use `--orthologs path/to/orthologs_1to1.csv`.

---

**Questions or tweaks?** Add labeling, facets, or alternative fits (e.g., robust regression) directly in the script as needed.

