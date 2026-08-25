#!/usr/bin/env python3
"""Normalize QC-local contig labels for PLINK and copy the matching FAI.

Port of the snpArcher qc module `prepare_plink_inputs` run block
(workflow/modules/qc/Snakefile): rewrites the reference FAI with the
QC-local plink_contig names and runs `bcftools annotate --rename-chrs`
on the pruned VCF.  Requires bcftools on PATH (the qc_subsample
environment).

Usage:
    prepare_plink_inputs.py <vcf> <contig_map.tsv> <ref.fai>
        <out.vcf> <out.fai> <log>
"""

import os
import subprocess
import sys
import tempfile


def load_contig_map(map_file):
    rename_map = {}
    with open(map_file) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        required = {"original_contig", "plink_contig"}
        if not required.issubset(set(header)):
            raise ValueError(f"Contig map is missing required column(s): {sorted(required)}")
        orig_i = header.index("original_contig")
        plink_i = header.index("plink_contig")
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            rename_map[fields[orig_i]] = fields[plink_i]
    return rename_map


def main():
    if len(sys.argv) != 7:
        print(__doc__)
        sys.exit(2)
    vcf, map_file, fai, out_vcf, out_fai, log = sys.argv[1:]

    rename_map = load_contig_map(map_file)

    with open(fai) as src, open(out_fai, "w") as dst:
        for line in src:
            fields = line.rstrip("\n").split("\t")
            if fields:
                fields[0] = rename_map.get(fields[0], fields[0])
            dst.write("\t".join(fields) + "\n")

    with tempfile.NamedTemporaryFile(
        "w",
        delete=False,
        dir=os.getcwd(),
        suffix=".rename.tsv",
    ) as rename_file:
        for original, plink_name in rename_map.items():
            rename_file.write(f"{original}\t{plink_name}\n")
        rename_path = rename_file.name

    try:
        with open(log, "wb") as log_handle:
            subprocess.run(
                ["bcftools", "annotate", "--rename-chrs", rename_path, vcf,
                 "-O", "z", "-o", out_vcf],
                stdout=log_handle, stderr=log_handle, check=False,
            )
            subprocess.run(
                ["bcftools", "index", out_vcf],
                stdout=log_handle, stderr=log_handle, check=False,
            )
    finally:
        os.unlink(rename_path)


if __name__ == "__main__":
    main()
