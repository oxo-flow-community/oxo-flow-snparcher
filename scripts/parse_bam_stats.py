#!/usr/bin/env python3
"""Parse samtools coverage + flagstat output into a per-sample JSON.

Port of the snpArcher `parse_bam_stats` run block
(workflow/rules/qc_metrics.smk): weighted-mean depth from the coverage
table, read counts from the TSV flagstat output. v2.3.0 keys the flagstat
lines by their description column (samtools labels the rows; positional
access breaks across samtools releases and on unpaired BAMs), and reports
"properly paired %" as None (NA in the combined report) when the library is
unpaired.

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

    # Parse flagstat -O tsv: QC-passed count, QC-failed count, description.
    flagstat = {}
    with open(flagstat_file) as f:
        for line in f:
            passed, _failed, description = line.rstrip("\n").split("\t")
            flagstat[description] = passed

    def flagstat_percent(description):
        value = flagstat[description].rstrip("%")
        return None if value == "N/A" else float(value)

    # Count primary records only. Secondary and supplementary records are
    # extra alignments of reads already counted, and would inflate the read
    # total and skew the mapping and duplicate rates.
    total_reads = int(flagstat["primary"])
    num_dups = int(flagstat["primary duplicates"])
    num_mapped = int(flagstat["primary mapped"])
    pct_mapped = flagstat_percent("primary mapped %") if total_reads > 0 else 0
    # None (NA in the report) for an unpaired BAM, which has no pairs.
    pct_proper_paired = (
        flagstat_percent("properly paired %") if total_reads > 0 else 0
    )

    out = {
        "sample": sample,
        "total_reads": total_reads,
        "num_mapped": num_mapped,
        "percent_mapped": pct_mapped,
        "num_duplicates": num_dups,
        "percent_duplicates": (
            num_dups / total_reads * 100 if total_reads > 0 else 0
        ),
        "percent_properly_paired": pct_proper_paired,
        "mean_depth": mean_depth,
        "covered_bases": covered_bases,
    }

    with open(output_file, "w") as f:
        json.dump(out, f, indent=2)


if __name__ == "__main__":
    main(sys.argv)
