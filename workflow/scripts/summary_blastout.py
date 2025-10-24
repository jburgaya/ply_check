#!/usr/bin/env python

def get_options():
    import argparse

    description = "Summarize blast output into output df"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("--out-dir",
                        help="Blast output dir")
    parser.add_argument("--out-file",
                        help="Output file to store results")

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    import os
    import glob
    import pandas as pd

    cols = ["allele", "contig", "pident", "length", "qstart", "qend", "sstart", "send"]

    records = []

    for filepath in glob.glob(os.path.join(options.out_dir, "*.tsv")):
        sample = os.path.basename(filepath).replace(".tsv", "")
        df = pd.read_csv(filepath, sep="\t", names=cols)

        if df.empty:
            # Empty BLAST file → ply-NP
            records.append({
                "sample": sample,
                "assigned_allele": "ply-NP",
                "pident": "",
                "alignment_length": "",
                "contig": "",
                "qstart": "",
                "qend": "",
                "sstart": "",
                "send": ""
            })
            continue

        best_hit = df.sort_values("pident", ascending=False).iloc[0]

        records.append({
            "sample": sample,
            "assigned_allele": best_hit["allele"],
            "pident": best_hit["pident"],
            "alignment_length": best_hit["length"],
            "contig": best_hit["contig"],
            "qstart": best_hit["qstart"],
            "qend": best_hit["qend"],
            "sstart": best_hit["sstart"],
            "send": best_hit["send"]
        })

    outfile = options.out_file
    if not outfile.endswith(".tsv"):
        outfile += ".tsv"

    result_df = pd.DataFrame(records)
    result_df.to_csv(outfile, sep="\t", index=False)
