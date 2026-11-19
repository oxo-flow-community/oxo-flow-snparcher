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

    # v2.3.0: the plink labels are always the positional qcctgN names, not
    # only when every contig already looks numeric.  QC-internal labels only
    # (the contig map drives the --rename-chrs/PLINK pair), so unifying them
    # removes the mixed-label class of PLINK contig mismatches.
    plink_contigs = [f"qcctg{i}" for i in range(1, len(contigs) + 1)]

    with open(output_path, "w") as handle:
        handle.write("original_contig\tplink_contig\tadmixture_id\n")
        for i, (original, plink) in enumerate(zip(contigs, plink_contigs), start=1):
            handle.write(f"{original}\t{plink}\t{i}\n")


if __name__ == "__main__":
    main()
