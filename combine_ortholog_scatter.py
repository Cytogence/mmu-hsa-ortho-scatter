#!/usr/bin/env python3
"""
Combine human & mouse DE results via 1:1 orthologs, export paired table, and draw a scatter.

Inputs:
  - human CSV:  first column = gene symbol; must contain 'log2FoldChange' and 'padj'
  - mouse CSV:  first column = gene symbol; must contain 'log2FoldChange' and 'padj'

By default, fetches 1:1 orthologs (human->mouse) from Ensembl BioMart.
Optionally, pass a precomputed ortholog CSV with columns: gene_human,gene_mouse.

Outputs:
  - --out-csv: CSV with columns:
      gene_mouse,gene_human,log2FC_mouse,log2FC_human,padj_mouse,padj_human
  - --out-plot: PNG scatter (x=mouse log2FC, y=human log2FC)
"""

import argparse, os, sys, time, html, io
import pandas as pd
import numpy as np
import requests
import matplotlib.pyplot as plt

BIOMART_URL = "https://www.ensembl.org/biomart/martservice"
# Alternative BioMart URLs to try if main one fails
BIOMART_MIRRORS = [
    "https://www.ensembl.org/biomart/martservice",
    "https://uswest.ensembl.org/biomart/martservice",
    "https://useast.ensembl.org/biomart/martservice",
    "https://asia.ensembl.org/biomart/martservice"
]

def read_deg_csv(path, species_label):
    df = pd.read_csv(path)
    if df.shape[1] < 3:
        raise ValueError(f"{species_label}: expected at least 3 columns (gene, log2FoldChange, padj); got {df.shape[1]}")
    # First column is the gene symbol
    gene_col = df.columns[0]
    # Find required columns (case-sensitive per your spec)
    if "log2FoldChange" not in df.columns or "padj" not in df.columns:
        raise ValueError(f"{species_label}: CSV must contain columns 'log2FoldChange' and 'padj'. Found: {list(df.columns)}")

    out = df[[gene_col, "log2FoldChange", "padj"]].copy()
    out.columns = ["gene", "log2FC", "padj"]
    # Clean symbols
    out["gene"] = out["gene"].astype(str).str.strip()
    out = out.dropna(subset=["gene"])
    return out

def chunked(iterable, n):
    chunk = []
    for x in iterable:
        chunk.append(x)
        if len(chunk) == n:
            yield chunk
            chunk = []
    if chunk:
        yield chunk

def fetch_orthologs_human_to_mouse(human_symbols, max_retries=10):
    """
    Query Ensembl BioMart for human->mouse one2one orthologs.
    Returns DataFrame with columns: gene_human, gene_mouse
    """
    if not human_symbols:
        return pd.DataFrame(columns=["gene_human","gene_mouse"])

    attrs = [
        "external_gene_name",
        "mmusculus_homolog_associated_gene_name",
        "mmusculus_homolog_orthology_type",
    ]
    
    frames = []
    total_batches = len(list(chunked(sorted(set(human_symbols)), 300)))
    
    for batch_num, batch in enumerate(chunked(sorted(set(human_symbols)), 300), 1):
        print(f"[info] Processing batch {batch_num}/{total_batches} ({len(batch)} genes)...", file=sys.stderr)
        
        values_xml = "".join(f"<Value>{html.escape(s)}</Value>" for s in batch)
        xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Filter name="hgnc_symbol" excluded="0">{values_xml}</Filter>
    <Filter name="with_mmusculus_homolog" excluded="0"/>
    {"".join(f'<Attribute name="{a}"/>' for a in attrs)}
  </Dataset>
</Query>"""
        
        # Try multiple mirrors and retry logic
        batch_success = False
        for mirror in BIOMART_MIRRORS:
            if batch_success:
                break
                
            for attempt in range(max_retries):
                try:
                    print(f"[info] Trying {mirror} (attempt {attempt + 1})...", file=sys.stderr)
                    
                    r = requests.post(mirror, data={"query": xml}, timeout=120)
                    
                    if r.status_code == 200 and r.text and not any(error in r.text.lower() for error in 
                        ["query error", "error:", "exception", "can't locate", "biomart::exception"]):
                        
                        # Check if we got actual data
                        lines = r.text.strip().split('\n')
                        if len(lines) > 1:  # More than just header
                            tmp = pd.read_csv(io.StringIO(r.text), sep="\t")
                            if not tmp.empty:
                                frames.append(tmp)
                                batch_success = True
                                print(f"[info] Batch {batch_num} successful with {len(tmp)} results", file=sys.stderr)
                                break
                    
                    # If we get here, the response wasn't good
                    if "Query ERROR" in r.text or "Exception" in r.text:
                        print(f"[warn] BioMart server error: {r.text[:200]}", file=sys.stderr)
                    else:
                        print(f"[warn] Unexpected response (HTTP {r.status_code}): {r.text[:100]}", file=sys.stderr)
                        
                except requests.exceptions.RequestException as e:
                    print(f"[warn] Network error on attempt {attempt + 1}: {e}", file=sys.stderr)
                except Exception as e:
                    print(f"[warn] Unexpected error on attempt {attempt + 1}: {e}", file=sys.stderr)
                
                if attempt < max_retries - 1:
                    wait_time = min(2 ** attempt, 60)  # Exponential backoff up to 60 seconds
                    print(f"[info] Waiting {wait_time} seconds before retry...", file=sys.stderr)
                    time.sleep(wait_time)
        
        if not batch_success:
            print(f"[error] Failed to fetch batch {batch_num} after trying all mirrors", file=sys.stderr)
            # Continue with other batches rather than failing completely

    if not frames:
        print("[error] No batches were successful. BioMart may be down.", file=sys.stderr)
        print("[error] Please try again later or provide a pre-computed ortholog file using --orthologs", file=sys.stderr)
        return pd.DataFrame(columns=["gene_human","gene_mouse"])

    raw = pd.concat(frames, ignore_index=True)
    print(f"[info] Retrieved {len(raw)} ortholog records from BioMart", file=sys.stderr)
    
    # Ensure expected columns exist
    for c in ["external_gene_name", "mmusculus_homolog_associated_gene_name", "mmusculus_homolog_orthology_type"]:
        if c not in raw.columns:
            raw[c] = pd.NA

    df = raw.loc[
        raw["mmusculus_homolog_orthology_type"] == "ortholog_one2one",
        ["external_gene_name", "mmusculus_homolog_associated_gene_name"]
    ].copy()

    df.columns = ["gene_human", "gene_mouse"]
    df = df.dropna().drop_duplicates()
    print(f"[info] Found {len(df)} one-to-one orthologs", file=sys.stderr)
    
    # Enforce strict 1:1 uniqueness on symbols
    dup_h = df["gene_human"].duplicated(keep=False)
    dup_m = df["gene_mouse"].duplicated(keep=False)
    df = df[~dup_h & ~dup_m].reset_index(drop=True)
    print(f"[info] After removing duplicates: {len(df)} unique 1:1 orthologs", file=sys.stderr)
    
    return df

def identify_genes_to_label(sig, x, y, n_labels=8):
    """
    Identify interesting genes to label based on extreme values and distances
    """
    try:
        # Calculate distances from origin
        distances = np.sqrt(x**2 + y**2)
        
        # Get indices of most extreme points
        extreme_indices = set()
        
        # Top and bottom in each dimension
        n_per_direction = max(1, n_labels // 4)
        
        # Highest/lowest x values
        extreme_indices.update(np.argsort(x)[-n_per_direction:])  # Highest x
        extreme_indices.update(np.argsort(x)[:n_per_direction])   # Lowest x
        
        # Highest/lowest y values  
        extreme_indices.update(np.argsort(y)[-n_per_direction:])  # Highest y
        extreme_indices.update(np.argsort(y)[:n_per_direction])   # Lowest y
        
        # Also add points with largest distance from origin
        extreme_indices.update(np.argsort(distances)[-n_per_direction:])
        
        # Convert to list and limit to n_labels
        label_indices = list(extreme_indices)[:n_labels]
        
        return label_indices
        
    except Exception as e:
        print(f"[warn] Could not identify genes to label: {e}", file=sys.stderr)
        return None

def ensure_dir(path):
    d = os.path.dirname(os.path.abspath(path))
    if d and not os.path.exists(d):
        os.makedirs(d, exist_ok=True)

def main():
    ap = argparse.ArgumentParser(description="Pair mouse & human DE by 1:1 orthologs and plot scatter.")
    ap.add_argument("--human", required=True, help="Human CSV (col1=gene, plus log2FoldChange & padj)")
    ap.add_argument("--mouse", required=True, help="Mouse CSV (col1=gene, plus log2FoldChange & padj)")
    ap.add_argument("--orthologs", help="Optional CSV with columns gene_human,gene_mouse (skip online fetch)")
    ap.add_argument("--out-csv", default="paired_mouse_human.csv", help="Output paired CSV path")
    ap.add_argument("--out-plot", default="mouse_human_scatter.png", help="Output plot PNG path")
    ap.add_argument("--padj-threshold", type=float, default=0.05, help="Significance cutoff for plotting")
    ap.add_argument("--sig-mode", choices=["both","any","human","mouse"], default="both",
                    help="Which padj must pass to include in plot (default: both)")
    args = ap.parse_args()

    # Read DE tables
    human = read_deg_csv(args.human, "human")
    mouse = read_deg_csv(args.mouse, "mouse")
    print(f"[info] Loaded {len(human)} human genes, {len(mouse)} mouse genes", file=sys.stderr)

    # Build or read ortholog table
    if args.orthologs:
        print(f"[info] Using pre-computed ortholog file: {args.orthologs}", file=sys.stderr)
        ortho = pd.read_csv(args.orthologs)
        if not {"gene_human","gene_mouse"}.issubset(ortho.columns):
            raise ValueError("orthologs CSV must contain columns: gene_human,gene_mouse")
        ortho = ortho[["gene_human","gene_mouse"]].dropna().drop_duplicates()
    else:
        print("[info] Fetching 1:1 orthologs from Ensembl BioMart...", file=sys.stderr)
        ortho = fetch_orthologs_human_to_mouse(human_symbols=set(human["gene"]))

    if ortho.empty:
        raise RuntimeError("No ortholog pairs found. Check gene symbols or provide --orthologs.")

    print(f"[info] Using {len(ortho)} ortholog pairs for matching", file=sys.stderr)

    # Merge: keep genes present in both DE tables and in the ortholog map
    merged = (ortho
              .merge(human.rename(columns={"gene":"gene_human","log2FC":"log2FC_human","padj":"padj_human"}),
                     on="gene_human", how="inner")
              .merge(mouse.rename(columns={"gene":"gene_mouse","log2FC":"log2FC_mouse","padj":"padj_mouse"}),
                     on="gene_mouse", how="inner"))

    if merged.empty:
        raise RuntimeError("No overlapping orthologs between your human and mouse DE tables after merge.")

    print(f"[info] Successfully matched {len(merged)} ortholog pairs between datasets", file=sys.stderr)

    # Write the paired CSV in your requested column order
    out_cols = ["gene_mouse","gene_human","log2FC_mouse","log2FC_human","padj_mouse","padj_human"]
    ensure_dir(args.out_csv)
    merged[out_cols].to_csv(args.out_csv, index=False)
    print(f"[ok] Wrote paired table: {args.out_csv}  (n={len(merged)})", file=sys.stderr)

    # Filter for plotting
    sig = merged.copy()
    thr = args.padj_threshold
    if args.sig_mode == "both":
        sig = sig[(sig["padj_mouse"] < thr) & (sig["padj_human"] < thr)]
    elif args.sig_mode == "any":
        sig = sig[(sig["padj_mouse"] < thr) | (sig["padj_human"] < thr)]
    elif args.sig_mode == "human":
        sig = sig[sig["padj_human"] < thr]
    elif args.sig_mode == "mouse":
        sig = sig[sig["padj_mouse"] < thr]

    if sig.empty:
        print("[warn] No points passed significance filter for plotting; skipping plot.", file=sys.stderr)
        return

    print(f"[info] Plotting {len(sig)} significant genes", file=sys.stderr)
    
    x = sig["log2FC_mouse"].astype(float).to_numpy()
    y = sig["log2FC_human"].astype(float).to_numpy()

    concordant_mask = (x > 0) & (y > 0) | (x < 0) & (y < 0)  # Same direction
    discordant_mask = (x > 0) & (y < 0) | (x < 0) & (y > 0)  # Opposite direction

    # Calculate statistics
    if len(x) >= 2:
        r = np.corrcoef(x, y)[0, 1]
    else:
        r = np.nan

    # Fit regression line
    try:
        slope, intercept = np.polyfit(x, y, 1)
    except Exception:
        slope, intercept = np.nan, np.nan

    # Create enhanced plot
    ensure_dir(args.out_plot)
    fig, ax = plt.subplots(figsize=(8, 7), dpi=150)
    
    # Create scatter plot with better styling
    ax.scatter(x[concordant_mask], y[concordant_mask], s=20, alpha=0.7, 
           c='lightgreen', edgecolors='white', linewidth=0.5, label='Concordant')
    ax.scatter(x[discordant_mask], y[discordant_mask], s=20, alpha=0.7, 
            c='lightcoral', edgecolors='white', linewidth=0.5, label='Discordant')
    
    # Add reference lines
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.axvline(0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    
    # Add regression line if valid
    if not np.isnan(slope):
        x_range = np.linspace(x.min(), x.max(), 100)
        y_pred = slope * x_range + intercept
        ax.plot(x_range, y_pred, color='red', linewidth=2, alpha=0.8)
    
    # Set labels and title with proper spacing
    ax.set_xlabel('Mouse log₂ Fold Change', fontweight='bold', fontsize=11)
    ax.set_ylabel('Human log₂ Fold Change', fontweight='bold', fontsize=11)
    ax.set_title('Cross-Species Differential Expression Correlation', 
                fontweight='bold', fontsize=12, pad=20)
    
    # Add statistics box in upper left corner
    stats_text = []
    if not np.isnan(r):
        stats_text.append(f'Pearson r = {r:.3f}')
        stats_text.append(f'R² = {r**2:.3f}')
    if not np.isnan(slope):
        stats_text.append(f'Slope = {slope:.3f}')
    stats_text.append(f'n = {len(sig)} genes')
    
    stats_str = '\n'.join(stats_text)
    ax.text(0.05, 0.95, stats_str, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', 
            facecolor='white', alpha=0.8, edgecolor='gray'))
    
    # Add gene labels for extreme points (8 most interesting genes)
    label_genes = identify_genes_to_label(sig, x, y, n_labels=20)
    
    if label_genes is not None:
        for idx in label_genes:
            gene_human = sig.iloc[idx]['gene_human']
            gene_mouse = sig.iloc[idx]['gene_mouse']
            
            # Use human gene name as primary, mouse as secondary
            if pd.notna(gene_human) and str(gene_human).strip():
                label = str(gene_human)
            elif pd.notna(gene_mouse) and str(gene_mouse).strip():
                label = str(gene_mouse)
            else:
                continue
                
            ax.annotate(label, (x[idx], y[idx]), 
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, ha='left', va='bottom',
                       bbox=dict(boxstyle='round,pad=0.2', facecolor='yellow', alpha=0.7),
                       arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    
    # Add legend for gene label types
    from matplotlib.patches import Rectangle
    legend_elements = [
        Rectangle((0,0),1,1, facecolor='lightgreen', alpha=0.8, edgecolor='darkgreen', 
                 label='Concordant genes (similar response)'),
        Rectangle((0,0),1,1, facecolor='yellow', alpha=0.8, edgecolor='orange',
                 label='Extreme genes (strongest responses)')
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9, framealpha=0.9)
    
    # Improve grid and layout
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    plt.tight_layout()
    
    # Save plot
    plt.savefig(args.out_plot, bbox_inches='tight', facecolor='white', dpi=150)
    plt.close()
    
    print(f"[ok] Wrote enhanced scatter: {args.out_plot}  (n points={len(sig)})", file=sys.stderr)

if __name__ == "__main__":
    main()