
import csv
import os

from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord

Entrez.email = "sanaparvin345@gmail.com"  # <-- CHANGE THIS (NCBI requires it)

ACCESSION = "NM_000518.5"

# Known HBB variants, written as c. positions (1-based, A of ATG = 1).
# Verify each one on ClinVar and cite it in your README.
VARIANTS = {
    "HbS (sickle) c.20A>T": (20, "A", "T"),
    "HbC c.19G>A": (19, "G", "A"),
    "HbE c.79G>A": (79, "G", "A"),
}

os.makedirs("results", exist_ok=True)


def fetch_record(accession):
    handle = Entrez.efetch(db="nucleotide", id=accession,
                           rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record


def get_cds(record):
    for feature in record.features:
        if feature.type == "CDS":
            return feature.extract(record.seq)
    raise ValueError("No CDS feature found")


def apply_mutation(cds, pos, ref, alt):
    if cds[pos - 1] != ref:
        raise ValueError(f"Expected {ref} at c.{pos}, found {cds[pos - 1]}")
    mutant = MutableSeq(cds)
    mutant[pos - 1] = alt
    return Seq(mutant)


def classify(wt_prot, mut_prot):
    """Return effect type plus protein change in HGVS-like notation."""
    for i, (a, b) in enumerate(zip(wt_prot, mut_prot), start=1):
        if a != b:
            if b == "*":
                return "nonsense", f"p.{a}{i}*"
            return "missense", f"p.{a}{i}{b}"
    return "synonymous", "p.="


def make_aligner():
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    return aligner


def main():
    print(f"Fetching {ACCESSION} from NCBI...")
    record = fetch_record(ACCESSION)
    cds = get_cds(record)
    wt_prot = cds.translate(to_stop=True)
    print(f"CDS length: {len(cds)} bp | Protein length: {len(wt_prot)} aa")

    SeqIO.write(SeqRecord(cds, id=ACCESSION, description="HBB CDS wild-type"),
                "results/hbb_wt.fasta", "fasta")

    aligner = make_aligner()
    rows = []
    with open("results/alignments.txt", "w") as aln_file:
        for name, (pos, ref, alt) in VARIANTS.items():
            mut_cds = apply_mutation(cds, pos, ref, alt)
            mut_prot = mut_cds.translate(to_stop=True)
            effect, change = classify(wt_prot, mut_prot)

            alignment = aligner.align(wt_prot, mut_prot)[0]
            identity = sum(a == b for a, b in zip(wt_prot, mut_prot)) \
                / len(wt_prot) * 100
            codon = (pos - 1) // 3 + 1

            rows.append({
                "variant": name,
                "cds_position": pos,
                "codon_number": codon,
                "ref": ref,
                "alt": alt,
                "effect": effect,
                "protein_change": change,
                "protein_identity_pct": round(identity, 2),
                "alignment_score": alignment.score,
            })
            aln_file.write(f"=== {name} ({change}, {effect}) ===\n")
            aln_file.write(str(alignment) + "\n\n")
            print(f"{name}: {effect} {change} | identity {identity:.1f}%")

    with open("results/variants.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print("Done. See the results/ folder.")

main()