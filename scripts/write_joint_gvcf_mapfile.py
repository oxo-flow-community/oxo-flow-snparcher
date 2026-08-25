#!/usr/bin/env python3
"""Write the GenomicsDB sample-name map from per-sample gVCF paths.

GenomicsDBImport requires a two-column TSV: sample name and gVCF path.
Sample names are derived from the `results/gvcfs/<sample>.g.vcf.gz` path
convention used by this workflow (upstream derives them from the sample
sheet; the port's gVCF paths carry the sample id).

Usage:
    write_joint_gvcf_mapfile.py <output_mapfile> <gvcf...>
"""

import os
import sys


def sample_from_gvcf(path):
    base = os.path.basename(path)
    for suffix in (".g.vcf.gz.tbi", ".g.vcf.gz", ".g.vcf.idx", ".g.vcf"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    raise ValueError(f"Unrecognized gVCF path (cannot derive sample id): {path}")


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    output = sys.argv[1]
    gvcfs = sys.argv[2:]
    with open(output, "w") as handle:
        for gvcf in gvcfs:
            sample = sample_from_gvcf(gvcf)
            handle.write(f"{sample}\t{gvcf}\n")


if __name__ == "__main__":
    main()
