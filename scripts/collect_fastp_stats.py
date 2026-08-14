#!/usr/bin/env python3
"""Aggregate fastp per-unit JSON stats into a per-sample JSON.

Port of the snpArcher `collect_fastp_stats` run block
(workflow/rules/fastq.smk): sums total reads before/after filtering
across all input units of a sample.

Usage: collect_fastp_stats.py <input.json>... <output.json>
"""
import json
import sys


def main(argv):
    if len(argv) < 2:
        raise SystemExit(f"usage: {argv[0]} <input.json>... <output.json>")
    input_files = argv[1:-1]
    output_file = argv[-1]

    unfiltered = 0
    pass_filter = 0

    for fn in input_files:
        with open(fn) as f:
            data = json.load(f)
        unfiltered += data["summary"]["before_filtering"]["total_reads"]
        pass_filter += data["summary"]["after_filtering"]["total_reads"]

    out = {
        "summary": {
            "before_filtering": {"total_reads": unfiltered},
            "after_filtering": {"total_reads": pass_filter},
        }
    }

    with open(output_file, "w") as f:
        json.dump(out, f, indent=2)


if __name__ == "__main__":
    main(sys.argv)
