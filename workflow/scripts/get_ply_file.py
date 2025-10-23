input_file = "ply_alleles.txt"
output_file = "ply_alleles.fasta"

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    sequence = ""
    header = ""

    for line in infile:
        line = line.strip()

        if line.startswith(">"):
            # Write the previous entry if exists
            if header and sequence:
                outfile.write(f"{header}\n{sequence}\n")
                sequence = ""

            # Extract ply-XX from header
            if "ply-" in line:
                ply_tag = line.split("ply-")[-1].split()[0].replace(",", "")
                header = f">ply-{ply_tag}"
            else:
                header = ">unknown"

        else:
            sequence += line

    # Write the last record
    if header and sequence:
        outfile.write(f"{header}\n{sequence}\n")

print(f"FASTA cleaned and saved to '{output_file}'")
