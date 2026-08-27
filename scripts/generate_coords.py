#!/usr/bin/env python3
"""Generate the QC coordinates file from optional sample metadata (lat/long).

Faithful port of upstream snpArcher v2.2 `modules/qc/Snakefile`
`generate_coords_file`: reads the optional `sample_metadata` CSV (columns
`sample_id`, `long`, `lat`), drops duplicate sample ids and rows without
coordinates, and writes a tab-separated file with no header. With no
metadata configured (or without lat/long columns) the file is written
empty, matching upstream's placeholder branch — the dashboard Rmd renders
the terrain map panel only when the file has rows AND a Google API key is
configured (`qc_google_api_key`).

Usage:
    generate_coords.py <sample_metadata_csv|''> <output_coords_txt>
"""

import sys

import pandas as pd


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    metadata_path, output = sys.argv[1], sys.argv[2]

    if not metadata_path:
        open(output, "w").close()
        return

    meta = pd.read_csv(metadata_path)
    if not {"sample_id", "long", "lat"}.issubset(meta.columns):
        open(output, "w").close()
        return

    out = meta[["sample_id", "long", "lat"]].copy()
    out.drop_duplicates("sample_id", inplace=True)
    out.dropna(subset=["long", "lat"], inplace=True)
    out.to_csv(output, index=False, sep="\t", header=False)


if __name__ == "__main__":
    main()
