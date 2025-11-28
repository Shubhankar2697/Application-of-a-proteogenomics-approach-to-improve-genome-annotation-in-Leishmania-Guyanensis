import pandas as pd
from Bio import SeqIO
import re
from pathlib import Path

species_files = {
 
    "L. guyanensis": (
        "L. guyanensis peptides.xlsx",
        "L. guyanensis hybrid 8066 proteins_070925.fasta"
    )

}

def clean_peptide(peptide):
    return re.sub(r'[^A-Z]', '', peptide)


def load_proteins(fasta_path):
    return [str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")]

# Match peptide sequence against proteins
def match_peptides(df, seq_column, protein_seqs):
    matched_rows, unmatched_rows = [], []
    for _, row in df.iterrows():
        raw_seq = str(row[seq_column])
        peptide = clean_peptide(raw_seq)
        found = any(peptide in protein for protein in protein_seqs)
        if found:
            matched_rows.append(row)
        else:
            unmatched_rows.append(row)
    return pd.DataFrame(matched_rows), pd.DataFrame(unmatched_rows)

# Process each species (If multiple input file)
for species, (xlsx_file, fasta_file) in species_files.items():
    print(f"Processing {species}")

    if not Path(xlsx_file).exists() or not Path(fasta_file).exists():
        print(f"Missing files for {species}, skipping.")
        continue

    df = pd.read_excel(xlsx_file)
    seq_col = next((col for col in df.columns if col.strip().lower() == "sequence"), None)

    if not seq_col:
        print(f"'Sequence' column not found in {xlsx_file}")
        continue

    protein_seqs = load_proteins(fasta_file)
    matched_df, unmatched_df = match_peptides(df, seq_col, protein_seqs)

    output_file = f"{species} mapped peptides.xlsx"
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        matched_df.to_excel(writer, sheet_name="Matched", index=False)
        unmatched_df.to_excel(writer, sheet_name="Unmatched", index=False)

    print(f"Done: {output_file}")

print("\n All species processed.")

