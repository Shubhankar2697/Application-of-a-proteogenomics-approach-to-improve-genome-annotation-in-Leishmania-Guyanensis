##Author: Shubhankar A Pawar, email: shubhankarpawar3@gmail.com; Dr.Harsh Pawar, Faculty Scientist, Institute of Bioinformatics, Bangalore email: harsh@ibioinformaics.org
##Six frame translation of L.guyanensis geneome MHOM/BR/75/M4147.
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable

standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
genetic_code = standard_table.forward_table
stop_codons = set(standard_table.stop_codons)

MIN_AA = 10
START_CODON = "ATG"

def scan_frame(seq_to_scan, offset, strand, seq_len, seq_base, accession, all_orfs):
    first_orf_in_frame = True
    pos = offset
    while pos + 3 <= seq_len:
        codon = seq_to_scan[pos:pos+3]
        if "N" in codon:
            pos += 1
            continue

        if first_orf_in_frame and codon != START_CODON:
            pos += 3
            continue

        # Start ORF here
        first_aa = "M" if codon == START_CODON else genetic_code.get(codon, "X")
        aa_seq = first_aa
        cur = pos + 3
        stop_found = False
        end_nt = None

        while cur + 3 <= seq_len:
            next_codon = seq_to_scan[cur:cur+3]
            if "N" in next_codon:
                end_nt = cur
                break
            if next_codon in stop_codons:
                end_nt = cur + 3
                stop_found = True
                break
            aa_seq += genetic_code.get(next_codon, "X")
            cur += 3

        if stop_found and len(aa_seq) >= MIN_AA:
            if strand == "+":
                g_start = pos + 1
                g_end = end_nt
            else:
                g_start = seq_len - end_nt + 1
                g_end = seq_len - pos

            all_orfs.append({
                "seq_base": seq_base,
                "accession": accession,
                "strand": strand,
                "start": g_start,
                "end": g_end,
                "aa": aa_seq,
                "start_codon": codon,
                "first_aa": first_aa
            })
            pos = end_nt
            first_orf_in_frame = False
            continue

        # No valid ORF starting here → move to next codon
        pos += 3

def run_orf_finder(input_fa, output_fa, output_log):
    all_orfs = []
    for rec in SeqIO.parse(input_fa, "fasta"):
        accession = rec.id
        seq_base = accession.replace(".", "_").replace("|", "_")
        seq_str = str(rec.seq).upper()
        L = len(seq_str)

        # Forward frames +1, +2, +3
        for offset in range(3):
            scan_frame(seq_str, offset, "+", L, seq_base, accession, all_orfs)

        # Reverse frames -1, -2, -3
        rc_str = str(Seq(seq_str).reverse_complement())
        for offset in range(3):
            scan_frame(rc_str, offset, "-", L, seq_base, accession, all_orfs)

    # Sort: accession → strand (+ before -) → start position
    all_orfs.sort(key=lambda x: (x["accession"], 0 if x["strand"] == "+" else 1, x["start"]))

    # Write output
    counters = {}
    with open(output_fa, "w") as fa, open(output_log, "w") as logf:
        for orf in all_orfs:
            key = (orf["accession"], orf["strand"])
            counters[key] = counters.get(key, 0) + 1
            n = counters[key]
            header = f"{orf['seq_base']}_{orf['strand']}_{orf['start']}-{orf['end']}_ORF{n}"
            fa.write(f">{header}\n{orf['aa']}\n")
            logf.write(f"{header} {orf['accession']} {orf['strand']} {orf['start']} {orf['end']} {orf['start_codon']} {orf['first_aa']}\n")

# ---- RUN ----
run_orf_finder("Lguyanesis.fasta", "Final_ORFs.fasta", "Final_ORFs.log")
print("DONE. ORFs saved to Final_ORFs.fasta and Final_ORFs.log")
