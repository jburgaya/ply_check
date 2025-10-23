#!/usr/bin/env python

def get_options():
    import argparse

    description = "Correct ply allele numbers in FASTA headers based on amino acid mutation profiles"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("--aa-changes",
                        required=True,
                        help="TSV file containing amino acid changes (aa_changes.tsv)")
    parser.add_argument("--out",
                        required=True,
                        help="Output")

    return parser.parse_args()


if __name__ == "__main__":
  import os
  import re
  import pandas as pd
  import numpy as np

  options = get_options()

  # Load amino acid changes table
  aa_changes = pd.read_csv(options.aa_changes, sep="\t")

  # reorder columns based on
  def extract_position(col):
      match = re.search(r'pos_(\d+)', col)
      return int(match.group(1)) if match else float('inf')

  # Keep 'sample_id' first, sort other columns by numeric position
  cols = aa_changes.columns.tolist()
  sorted_cols = ['sample_id'] + sorted([c for c in cols if c != 'sample_id'], key=extract_position)

  # Reorder DataFrame
  aa_changes = aa_changes[sorted_cols]
  aa_changes.head(2)

  # extract sampleid and ply columns
  # Extract sampleid (digits before _ply-) or full sample_id if format is different
  aa_changes["sampleid"] = aa_changes["sample_id"].str.extract(r'^(\d+)_ply-\d+')
  aa_changes["sampleid"] = aa_changes["sampleid"].fillna(aa_changes["sample_id"])  # fallback if no match

  # Extract ply allele (e.g. ply-5)
  aa_changes["ply_allele_nucl"] = aa_changes["sample_id"].str.extract(r'(ply-\d+)')

  # Reorder columns: sample_id, sampleid, ply_allele, then the rest
  first_cols = ["sample_id", "sampleid", "ply_allele_nucl"]
  remaining_cols = [col for col in aa_changes.columns if col not in first_cols]
  aa_changes = aa_changes[first_cols + remaining_cols]

  aa_changes.insert(aa_changes.columns.get_loc("ply_allele_nucl") + 1, "ply_allele_aa", "")

  # Identify amino acid mutation columns
  aa_cols = aa_changes.columns[aa_changes.columns.get_loc("ply_allele_aa") + 1:]

  # Extract reference rows (e.g., "ply-2", "ply-3", ...)
  ref_df  = aa_changes[aa_changes["sample_id"].str.match(r"^ply-\d+$")].copy()

  # Target only non-reference rows
  sample_df = aa_changes[~aa_changes["sample_id"].str.match(r"^ply-\d+$")].copy()

  for idx, sample_row in sample_df.iterrows():
      sample_aa = sample_row[aa_cols]

      # Loop through reference rows
      matched = False
      for _, ref_row in ref_df.iterrows():
          ref_aa = ref_row[aa_cols]

          # Only compare positions where the sample has non-NaN values
          mask = sample_aa.notna()

          # Enforce: all values must match exactly in those positions
          if mask.any() and sample_aa[mask].equals(ref_aa[mask]):
              # Set ply_allele_aa to the matching reference allele (ply_allele_nucl)
              aa_changes.at[idx, "ply_allele_aa"] = ref_row["ply_allele_nucl"]
              matched = True
              break

      if not matched:
          aa_changes.at[idx, "ply_allele_aa"] = ""  # or np.nan
  # save aa changes
  aa_changes.to_csv(options.out, sep="\t", index=None)
