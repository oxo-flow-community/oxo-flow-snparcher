#!/usr/bin/env python3
"""Write the GenomicsDB sample-name map from per-sample gVCF paths.

GenomicsDBImport requires a two-column TSV: sample name and gVCF path.
Sample names are derived from the `results/gvcfs/.../<sample>.g.vcf[.gz]`
path convention used by this workflow (upstream derives them from the sample
sheet; the port's gVCF paths carry the sample id).

The long-contig chain hands the rule its ordering-anchor inputs as well (the
durable .gz + CSI archive, the GATK .idx companion), so the input list can
carry several files per sample. One line is written per sample, picking the
path GATK must read: the plain work gVCF when present (long mode keeps GATK
on uncompressed input; upstream get_gatk_joint_gvcf), else the archived .gz
(short mode), never an index file.

Usage:
    write_joint_gvcf_mapfile.py <output_mapfile> <gvcf-or-companion...>
"""

import os
import sys

# Preference order for the representative path of one sample (v2.3.0
# get_gatk_joint_gvcf: GATK reads the plain work gVCF in long-contig mode and
# the archive otherwise). Index companions are never emitted.
PREFERRED_SUFFIXES = (".g.vcf", ".g.vcf.gz", ".g.vcf.gz.tbi", ".g.vcf.gz.csi")
IGNORED_SUFFIXES = (".g.vcf.idx",)


def sample_from_gvcf(path):
    base = os.path.basename(path)
    for suffix in (".g.vcf.gz.tbi", ".g.vcf.gz.csi", ".g.vcf.gz.idx",
                   ".g.vcf.idx", ".g.vcf.gz", ".g.vcf"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    raise ValueError(f"Unrecognized gVCF path (cannot derive sample id): {path}")


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    output = sys.argv[1]
    gvcfs = sys.argv[2:]
    # sample -> ranked candidate list, in first-seen order of samples.
    order = []
    candidates = {}
    for gvcf in gvcfs:
        base = os.path.basename(gvcf)
        if base.endswith(IGNORED_SUFFIXES):
            continue
        sample = sample_from_gvcf(gvcf)
        if sample not in candidates:
            candidates[sample] = []
            order.append(sample)
        candidates[sample].append(gvcf)
    with open(output, "w") as handle:
        for sample in order:
            paths = candidates[sample]

            def rank(path):
                base = os.path.basename(path)
                for idx, suffix in enumerate(PREFERRED_SUFFIXES):
                    if base.endswith(suffix):
                        return idx
                return len(PREFERRED_SUFFIXES)

            handle.write(f"{sample}\t{min(paths, key=rank)}\n")


if __name__ == "__main__":
    main()
