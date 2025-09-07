#!/usr/bin/env python3
"""
Download and prepare human-mouse ortholog file from Allen Institute
"""
import pandas as pd
import requests
import sys

def download_and_prepare_orthologs():
    """Download Allen Institute ortholog file and extract human-mouse pairs"""
    
    url = "https://github.com/AllenInstitute/GeneOrthology/raw/main/csv/mouse_human_marmoset_macaque_orthologs_20231113.csv"
    
    print("Downloading ortholog file from Allen Institute...", file=sys.stderr)
    
    try:
        # Download the file
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        
        # Read into pandas
        from io import StringIO
        df = pd.read_csv(StringIO(response.text))
        
        print(f"Downloaded {len(df)} ortholog records", file=sys.stderr)
        
        # Extract human and mouse columns
        # The Allen Institute file has columns: human_Symbol, mouse_Symbol
        if 'human_Symbol' not in df.columns or 'mouse_Symbol' not in df.columns:
            print(f"Expected columns not found. Available columns: {list(df.columns)}", file=sys.stderr)
            return None
            
        # Create clean ortholog mapping
        ortho = df[['human_Symbol', 'mouse_Symbol']].copy()
        ortho.columns = ['gene_human', 'gene_mouse']  # Rename to match your script's expected format
        
        # Remove rows with missing values
        ortho = ortho.dropna()
        
        # Remove duplicates
        ortho = ortho.drop_duplicates()
        
        print(f"Extracted {len(ortho)} clean human-mouse ortholog pairs", file=sys.stderr)
        
        # Save to file
        output_file = "human_mouse_orthologs.csv"
        ortho.to_csv(output_file, index=False)
        print(f"Saved ortholog file as: {output_file}", file=sys.stderr)
        
        return output_file
        
    except Exception as e:
        print(f"Error downloading orthologs: {e}", file=sys.stderr)
        return None

if __name__ == "__main__":
    output_file = download_and_prepare_orthologs()
    if output_file:
        print(f"\nSuccess! Use this file with your script:")
        print(f"--orthologs {output_file}")
    else:
        print("Failed to download orthologs. Try manual download.")