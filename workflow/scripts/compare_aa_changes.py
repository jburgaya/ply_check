#!/usr/bin/env python

import argparse
from Bio import SeqIO
import pandas as pd

def get_options():
    parser = argparse.ArgumentParser(description="Compare AA sequences to reference and list changes")

    parser.add_argument("--aln", required=True, help="Amino acid alignment FASTA file")
    parser.add_argument("--ref", default="ply-1", help="Reference sequence ID (default: ply-1)")
    parser.add_argument("--output", required=True, help="Output file")

    return parser.parse_args()

def main():
    options = get_options()

    records = list(SeqIO.parse(options.aln, "fasta"))
    ref_record = next((r for r in records if r.id == options.ref), None)
    if ref_record is None:
        raise ValueError(f"Reference ID '{options.ref}' not found in alignment.")

    ref_seq = str(ref_record.seq)
    alignment_length = len(ref_seq)

    changes = []

    for record in records:
        if record.id == options.ref:
            continue
        seq = str(record.seq)
        sample_changes = {"sample_id": record.id}
        for i in range(alignment_length):
            ref_aa = ref_seq[i]
            alt_aa = seq[i]
            if alt_aa != ref_aa:
                if ref_aa == "-":
                    col_name = f"pos_{i+1}_ins"
                elif alt_aa == "-":
                    col_name = f"pos_{i+1}_{ref_aa}_del"
                else:
                    col_name = f"pos_{i+1}_{ref_aa}"
                sample_changes[col_name] = alt_aa
        changes.append(sample_changes)

    df = pd.DataFrame(changes).fillna("")
    df.to_csv(options.output, index=False, sep="\t")
    print(f"Wrote AA change matrix to {options.output}")

if __name__ == "__main__":
    main()
