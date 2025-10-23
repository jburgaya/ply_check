#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_options():
    parser = argparse.ArgumentParser(description="Update FASTA headers based on mapping file (4th column)")
    parser.add_argument("--input-fasta", required=True, help="Input FASTA file")
    parser.add_argument("--mapping-file", required=True, help="Mapping file (tab-delimited, 4th column used)")
    parser.add_argument("--output-fasta", required=True, help="Output FASTA file with updated headers")
    return parser.parse_args()

def load_mapping(mapping_file):
    mapping = {}
    with open(mapping_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                prefix = parts[1]  # the column that matches the prefix before ply
                new_ply = parts[3]  # the new ply value
                mapping[prefix] = new_ply
    return mapping

def update_fasta_headers(input_fasta, mapping, output_fasta):
    updated_records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        if "_" in record.id and "ply-" in record.id.split("_")[1]:
            prefix, ply_part = record.id.split("_", 1)
            if prefix in mapping:
                record.id = f"{prefix}_{mapping[prefix]}"
                record.description = ""
        updated_records.append(record)
    SeqIO.write(updated_records, output_fasta, "fasta")

if __name__ == "__main__":
    options = get_options()
    mapping = load_mapping(options.mapping_file)
    update_fasta_headers(options.input_fasta, mapping, options.output_fasta)
