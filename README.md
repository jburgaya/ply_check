# plyCheck

Annotate ply alleles

``db/ply_alleles.txt``: allele sequences extracted from NCBI

``db/ply_alleles.fasta``: allele sequences to be used in the pipeline

# Usage

Prepare input files in ``config/config.yml`` file.

Create output directories ``òut/blast`` and ``out/logs``

Then  run snakemake as:

``snakemake --use-conda --cores 4 -R all -p``

# Outputs

| output  | description |
|----------|----------|
| ``blastoutput.tsv``    | Merged output from blast with statistics and initial ply allele assignation    |
| ``alleles.fna``    | File with extracted ply genes    |
| ``aa_unaln.fasta``    | File with extracted ply aa and reference aa   |
| ``aa_aln.fasta``    | Aligned  ``aa_unaln.fasta``    |
| ``aa_aln_corrected.fasta``    | Corrected ply allele assignment    |
| ``aa_changes.tsv``    | Initial ply allele assignment    |
| ``aa_changes_corrected.tsv``    | Corrected ply allele assignment    |

