#!/usr/bin/env python3
"""Combine per-sample QC stats into the final TSV report.

Port of the snpArcher `combine_qc_metrics` run block
(workflow/rules/qc_metrics.smk). Input files are per-sample fastp and
BAM stat JSONs (distinguished by path); output is the cohort
results/qc_metrics/qc_report.tsv.

Usage: combine_qc_metrics.py <input.json>... <output.tsv>
"""
import json
import os
import sys


def main(argv):
    if len(argv) < 2:
        raise SystemExit(
            f"usage: {argv[0]} <input.json>... <output.tsv>"
        )
    input_files = argv[1:-1]
    output_file = argv[-1]

    # Load fastp stats
    fastp_stats = {}
    # Load bam stats
    bam_stats = {}

    for fn in input_files:
        # gvcf cohorts produce neither fastp nor bam metrics — aggregate
        # whatever exists (upstream's qc module iterates existing files).
        if not os.path.exists(fn):
            continue
        sample = os.path.basename(fn).replace(".json", "")
        with open(fn) as f:
            data = json.load(f)
        if "/fastp/" in fn.replace("\\", "/"):
            fastp_stats[sample] = data
        else:
            bam_stats[sample] = data

    # Build header
    header = [
        "sample",
        "reads_before_filtering",
        "reads_after_filtering",
        "fraction_passed",
        "total_reads",
        "percent_mapped",
        "num_duplicates",
        "percent_duplicates",
        "percent_properly_paired",
        "mean_depth",
        "covered_bases",
    ]

    # Write report
    with open(output_file, "w") as f:
        f.write("\t".join(header) + "\n")

        for sample in bam_stats:
            row = [sample]

            # Fastp stats (may not exist for bam input type)
            if sample in fastp_stats:
                fs = fastp_stats[sample]
                row.extend(
                    [
                        str(fs["summary"]["before_filtering"]["total_reads"]),
                        str(fs["summary"]["after_filtering"]["total_reads"]),
                        f"{fs['summary']['after_filtering']['total_reads'] / fs['summary']['before_filtering']['total_reads']:.4f}",
                    ]
                )
            else:
                row.extend(["NA", "NA", "NA"])

            # BAM stats
            bam = bam_stats[sample]
            row.extend(
                [
                    str(bam["total_reads"]),
                    f"{bam['percent_mapped']:.2f}",
                    str(bam["num_duplicates"]),
                    f"{bam['percent_duplicates']:.2f}",
                    f"{bam['percent_properly_paired']:.2f}",
                    f"{bam['mean_depth']:.2f}",
                    str(bam["covered_bases"]),
                ]
            )

            f.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main(sys.argv)
