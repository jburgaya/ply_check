#!/usr/bin/env python

import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess

def get_options():
    parser = argparse.ArgumentParser(description="Translate and align allele sequences")
    parser.add_argument("--db-fasta", required=True, help="Path to reference allele FASTA (nucleotide)")
    parser.add_argument("--assigned-fasta", required=True, help="Path to assigned allele FASTA (nucleotide)")
    parser.add_argument("--output-aa-fasta", required=True, help="Output combined amino acid FASTA before alignment")
    parser.add_argument("--output-aln", required=True, help="Path to save aligned protein sequences (FASTA format)")
    return parser.parse_args()

def translate_records(nuc_fasta):
    aa_records = []
    for record in SeqIO.parse(nuc_fasta, "fasta"):
        aa_seq = record.seq.translate(to_stop=True)
        aa_record = SeqRecord(aa_seq, id=record.id, description="")
        aa_records.append(aa_record)
    return aa_records

if __name__ == "__main__":
    options = get_options()

    # Translate both fasta files
    db_aa = translate_records(options.db_fasta)
    assigned_aa = translate_records(options.assigned_fasta)

    all_aa = db_aa + assigned_aa

    # Write unaligned protein sequences
    SeqIO.write(all_aa, options.output_aa_fasta, "fasta")

    # Align using MAFFT (must be installed or in environment)
    subprocess.run([
        "mafft",
        "--auto",
        options.output_aa_fasta
    ], stdout=open(options.output_aln, "w"))
