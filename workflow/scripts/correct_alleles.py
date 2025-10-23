#!/usr/bin/env python

def get_options():
    import argparse

    description = "Correct ply allele numbers in FASTA headers based on amino acid mutation profiles"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("--aa-changes",
                        required=True,
                        help="TSV file containing amino acid changes (aa_changes.tsv)")
    parser.add_argument("--input-fasta",
                        required=True,
                        help="FASTA file with sequences to correct (headers like sampleid_ply-X)")
    parser.add_argument("--output-fasta",
                        required=True,
                        help="Output FASTA file with corrected headers")
    parser.add_argument("--summary-out",
                        required=False,
                        help="Optional output CSV with old/new allele mapping")
    parser.add_argument("--update-aa-changes",
                        action="store_true",
                        help="If set, writes an updated aa_changes file including the ply_allele_aa column")

    return parser.parse_args()


if __name__ == "__main__":
    import os
    import re
    import pandas as pd
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    options = get_options()

    # === Check files exist ===
    if not os.path.exists(options.aa_changes):
        raise FileNotFoundError(f"❌ Amino acid changes file not found: {options.aa_changes}")
    if not os.path.exists(options.input_fasta):
        raise FileNotFoundError(f"❌ Input FASTA file not found: {options.input_fasta}")

    print(f"🔹 Reading AA changes from: {options.aa_changes}")
    print(f"🔹 Reading FASTA from: {options.input_fasta}")

    # === Load amino acid change table ===
    aa_changes = pd.read_csv(options.aa_changes, sep="\t")

    # --- Reorder columns numerically by position ---
    def extract_position(col):
        match = re.search(r'pos_(\d+)', col)
        return int(match.group(1)) if match else float('inf')

    cols = aa_changes.columns.tolist()
    sorted_cols = ['sample_id'] + sorted(
        [c for c in cols if c != 'sample_id'],
        key=extract_position
    )
    aa_changes = aa_changes[sorted_cols]

    # --- Extract sampleid and allele ---
    aa_changes["sampleid"] = aa_changes["sample_id"].str.extract(r'^(\d+)_ply-\d+')
    aa_changes["sampleid"] = aa_changes["sampleid"].fillna(aa_changes["sample_id"])
    aa_changes["ply_allele_nucl"] = aa_changes["sample_id"].str.extract(r'(ply-\d+)')

    # Insert placeholder for corrected AA allele
    aa_changes.insert(aa_changes.columns.get_loc("ply_allele_nucl") + 1, "ply_allele_aa", "")

    # --- Identify amino acid columns ---
    aa_cols = aa_changes.columns[aa_changes.columns.get_loc("ply_allele_aa") + 1:]

    # --- Separate reference and sample rows ---
    ref_df = aa_changes[aa_changes["sample_id"].str.match(r"^ply-\d+$")].copy()
    sample_df = aa_changes[~aa_changes["sample_id"].str.match(r"^ply-\d+$")].copy()

    # === Match sample AA profiles to reference alleles ===
    print("🔹 Matching sample amino acid profiles to reference alleles...")
    for idx, sample_row in sample_df.iterrows():
        sample_aa = sample_row[aa_cols]
        matched = False

        for _, ref_row in ref_df.iterrows():
            ref_aa = ref_row[aa_cols]
            mask = sample_aa.notna()
            if mask.any() and sample_aa[mask].equals(ref_aa[mask]):
                aa_changes.at[idx, "ply_allele_aa"] = ref_row["ply_allele_nucl"]
                matched = True
                break

        if not matched:
            aa_changes.at[idx, "ply_allele_aa"] = ""  # no exact match

    # Merge sample + reference back
    aa_changes = pd.concat([ref_df, aa_changes.loc[sample_df.index]])

    # === Build header mapping ===
    mapping = {}
    for _, row in aa_changes.iterrows():
        if re.match(r'^\d+_ply-\d+$', row["sample_id"]):
            old_header = row["sample_id"]
            new_allele = row["ply_allele_aa"] if row["ply_allele_aa"] else row["ply_allele_nucl"]
            new_header = f"{row['sampleid']}_{new_allele}"
            mapping[old_header] = new_header

    # === Rename headers in FASTA ===
    print(f"🔹 Writing corrected FASTA to: {options.output_fasta}")
    updated = 0
    with open(options.output_fasta, "w") as out_f:
        for record in SeqIO.parse(options.input_fasta, "fasta"):
            old_id = record.id
            if old_id in mapping and mapping[old_id] != old_id:
                record.id = mapping[old_id]
                record.description = ""
                updated += 1
            SeqIO.write(record, out_f, "fasta")

    print(f"✅ FASTA correction complete — {updated} headers updated")

    # === Optional: Write updated AA changes file ===
    if options.update_aa_changes:
        aa_changes_out = os.path.splitext(options.aa_changes)[0] + "_updated.tsv"
        aa_changes.to_csv(aa_changes_out, sep="\t", index=False)
        print(f"✅ Updated amino acid change file written to: {aa_changes_out}")

    # === Optional: Write mapping summary ===
    if options.summary_out:
        mapping_df = pd.DataFrame([
            {"old_header": k, "new_header": v, "changed": k != v}
            for k, v in mapping.items()
        ])
        mapping_df.to_csv(options.summary_out, index=False)
        print(f"✅ Header mapping summary written to: {options.summary_out}")
