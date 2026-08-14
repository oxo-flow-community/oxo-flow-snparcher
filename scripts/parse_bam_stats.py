#!/usr/bin/env python3
"""Parse samtools coverage + flagstat output into a per-sample JSON.

Port of the snpArcher `parse_bam_stats` run block
(workflow/rules/qc_metrics.smk): weighted-mean depth from the coverage
table, read counts from the TSV flagstat output.

Usage: parse_bam_stats.py <coverage.txt> <flagstat.txt> <output.json>
"""
import json
import sys


def main(argv):
    if len(argv) != 4:
        raise SystemExit(
            f"usage: {argv[0]} <coverage.txt> <flagstat.txt> <output.json>"
        )
    coverage_file, flagstat_file, output_file = argv[1:]

    # Sample name derived from the coverage filename ({sample}_coverage.txt)
    sample = coverage_file.replace("\\", "/").split("/")[-1]
    if sample.endswith("_coverage.txt"):
        sample = sample[: -len("_coverage.txt")]

    # Parse coverage - weighted average across scaffolds
    num_sites = []
    depths = []
    covered_bases = 0

    with open(coverage_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            start, end = int(fields[1]), int(fields[2])
            num_sites.append(end - start + 1)
            depths.append(float(fields[6]))
            covered_bases += int(fields[4])

    total_sites = sum(num_sites)
    mean_depth = (
        sum(d * n / total_sites for d, n in zip(depths, num_sites))
        if total_sites > 0
        else 0
    )

    # Parse flagstat
    with open(flagstat_file) as f:
        lines = f.readlines()

    total_reads = int(lines[0].split()[0])
    num_dups = int(lines[4].split()[0])
    num_mapped = int(lines[6].split()[0])
    pct_mapped = (
        float(lines[7].split()[0].strip("%")) if total_reads > 0 else 0
    )
    pct_proper_paired = (
        float(lines[14].split()[0].strip("%")) if total_reads > 0 else 0
    )

    out = {
        "sample": sample,
        "total_reads": total_reads,
        "num_mapped": num_mapped,
        "percent_mapped": pct_mapped,
        "num_duplicates": num_dups,
        "percent_duplicates": (num_dups / total_reads * 100)
        if total_reads > 0
        else 0,
        "percent_properly_paired": pct_proper_paired,
        "mean_depth": mean_depth,
        "covered_bases": covered_bases,
    }

    with open(output_file, "w") as f:
        json.dump(out, f, indent=2)


if __name__ == "__main__":
    main(sys.argv)
