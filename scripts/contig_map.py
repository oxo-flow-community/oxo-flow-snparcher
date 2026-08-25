#!/usr/bin/env python3
"""Generate the QC-local contig map for PLINK and ADMIXTURE.

Port of the snpArcher qc module `contig_map` run block
(workflow/modules/qc/Snakefile).  The upstream block uses pandas; this
port keeps the same logic with the standard library.

Usage: contig_map.py <ref.fai> <contig_map.tsv>
"""

import sys


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    fai_path, output_path = sys.argv[1], sys.argv[2]

    with open(fai_path) as handle:
        contigs = [line.split("\t", 1)[0] for line in handle if line.strip()]
    if not contigs:
        raise ValueError(f"Reference FAI is empty: {fai_path}")

    all_numeric = all(c.replace(".", "", 1).isdigit() for c in contigs)
    if all_numeric:
        plink_contigs = [f"qcctg{i}" for i in range(1, len(contigs) + 1)]
    else:
        plink_contigs = list(contigs)

    with open(output_path, "w") as handle:
        handle.write("original_contig\tplink_contig\tadmixture_id\n")
        for i, (original, plink) in enumerate(zip(contigs, plink_contigs), start=1):
            handle.write(f"{original}\t{plink}\t{i}\n")


if __name__ == "__main__":
    main()
