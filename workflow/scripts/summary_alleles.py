#!/usr/bin/env python

def get_options():
    import argparse

    description = "Extract assigned allele sequences from genomes based on BLAST hits"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("--hits-file",
                        required=True,
                        help="Output file from summary_blast")
    parser.add_argument("--fasta-dir",
                        required=True,
                        help="Directory containing genome fasta files")
    parser.add_argument("--output-fasta",
                        required=True,
                        help="Output fasta file to save extracted allele sequences")

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    import os
    import pandas as pd
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    hits_df = pd.read_csv(options.hits_file, sep="\t")

    records = []
    skipped = 0

    for _, row in hits_df.iterrows():
        # Check for missing data
        if pd.isna(row["sstart"]) or pd.isna(row["send"]) or pd.isna(row["contig"]):
            print(f"⚠️ Skipping sample {row.get('sample', 'UNKNOWN')} due to missing BLAST fields.")
            skipped += 1
            continue

        sample = str(row["sample"])
        contig = row["contig"]
        start = int(row["sstart"])
        end = int(row["send"])
        allele = row["assigned_allele"]

        fasta_path = os.path.join(options.fasta_dir, sample + ".fasta")

        if not os.path.isfile(fasta_path):
            print(f"⚠️ Skipping sample {sample}: FASTA file not found at {fasta_path}")
            skipped += 1
            continue

        contigs = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

        if contig not in contigs:
            print(f"⚠️ Skipping sample {sample}: Contig '{contig}' not found in {fasta_path}")
            skipped += 1
            continue

        seq = contigs[contig].seq
        if start < end:
            allele_seq = seq[start-1:end]  # Convert 1-based to 0-based indexing
        else:
            allele_seq = seq[end-1:start].reverse_complement()

        seq_id = f"{sample}_{allele}"
        records.append(SeqRecord(allele_seq, id=seq_id, description=""))

    if records:
        SeqIO.write(records, options.output_fasta, "fasta")
        print(f"✅ Extracted {len(records)} allele sequences to {options.output_fasta}")
    else:
        print("❌ No sequences extracted.")

    if skipped > 0:
        print(f"⚠️ Skipped {skipped} entries due to errors or missing data.")
