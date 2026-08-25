#!/usr/bin/env python3
"""Write one sample id per line from a comma-separated sample list.

Port of the postprocess `filter_individuals` rule: upstream writes the
samples that are not excluded in the metadata sheet; the port's sample
sheet is the [[sample_groups]] list, and this script excludes nothing
(metadata-driven exclusion is not ported).

Usage:
    write_include_samples.py <output_file> <comma,separated,samples>
"""

import sys


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    output = sys.argv[1]
    samples = sys.argv[2]
    with open(output, "w") as handle:
        for sample in samples.split(","):
            sample = sample.strip()
            if sample:
                handle.write(sample + "\n")


if __name__ == "__main__":
    main()
