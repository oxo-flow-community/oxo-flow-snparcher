#!/usr/bin/env python3
"""Compute per-individual depth, missingness, het, and filter summary.

Port of the snpArcher qc module `vcftools_individuals` run block
(workflow/modules/qc/Snakefile): runs vcftools FILTER-summary/het/
missing-indv, computes per-sample SNP depth from FORMAT/DP (vcftools
--depth) or FORMAT/AD (summed AD alleles) or writes NA, then writes the
min_depth sample-include list.  Requires vcftools on PATH (the
qc_vcftools environment).

Usage:
    vcftools_individuals.py <vcf> <prefix> <min_depth> <log>
        <out.idepth> <out.imiss> <out.samps.txt> <out.summ> <out.het>
"""

import gzip
import math
import os
import subprocess
import sys


def _open_vcf_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _read_vcf_header(path):
    format_ids = set()
    samples = []
    with _open_vcf_text(path) as handle:
        for line in handle:
            if line.startswith("##FORMAT=<ID="):
                format_id = line.split("##FORMAT=<ID=", 1)[1].split(",", 1)[0]
                format_ids.add(format_id)
            elif line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                break
    if not samples:
        raise ValueError(f"Missing #CHROM header or sample columns in VCF: {path}")
    return format_ids, samples


def _finite_float(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def _normalize_depth_file(path):
    rows = []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split()
        required = {"INDV", "N_SITES", "MEAN_DEPTH"}
        if not required.issubset(set(header)):
            raise ValueError(f"Depth file is missing required columns: {path}")
        indv_i = header.index("INDV")
        nsites_i = header.index("N_SITES")
        mean_i = header.index("MEAN_DEPTH")
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split()
            mean_depth = _finite_float(fields[mean_i]) if len(fields) > mean_i else None
            rows.append((
                fields[indv_i],
                fields[nsites_i],
                "NA" if mean_depth is None else f"{mean_depth:.6g}",
            ))
    with open(path, "w") as handle:
        handle.write("INDV\tN_SITES\tMEAN_DEPTH\n")
        for indv, nsites, mean_depth in rows:
            handle.write(f"{indv}\t{nsites}\t{mean_depth}\n")


def _ad_depth(ad_value):
    if ad_value in {"", "."}:
        return None
    values = ad_value.split(",")
    if not values:
        return None
    depths = []
    for value in values:
        depth = _finite_float(value)
        if depth is None:
            return None
        depths.append(depth)
    return sum(depths)


def _write_ad_depth(vcf_path, depth_path, samples):
    sums = {sample: 0.0 for sample in samples}
    counts = {sample: 0 for sample in samples}
    with _open_vcf_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            format_fields = fields[8].split(":")
            if "AD" not in format_fields:
                continue
            ad_i = format_fields.index("AD")
            for sample, genotype in zip(samples, fields[9:]):
                genotype_fields = genotype.split(":")
                if len(genotype_fields) <= ad_i:
                    continue
                depth = _ad_depth(genotype_fields[ad_i])
                if depth is None:
                    continue
                sums[sample] += depth
                counts[sample] += 1
    with open(depth_path, "w") as handle:
        handle.write("INDV\tN_SITES\tMEAN_DEPTH\n")
        for sample in samples:
            if counts[sample] == 0:
                handle.write(f"{sample}\t0\tNA\n")
            else:
                handle.write(f"{sample}\t{counts[sample]}\t{sums[sample] / counts[sample]:.6g}\n")


def _write_empty_depth(depth_path, samples):
    with open(depth_path, "w") as handle:
        handle.write("INDV\tN_SITES\tMEAN_DEPTH\n")
        for sample in samples:
            handle.write(f"{sample}\t0\tNA\n")


def _write_sample_filter(depth_path, samples, samps_path, min_depth, log_path):
    passing = []
    finite_depths = 0
    with open(depth_path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        indv_i = header.index("INDV")
        mean_i = header.index("MEAN_DEPTH")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= mean_i:
                continue
            depth = _finite_float(fields[mean_i])
            if depth is None:
                continue
            finite_depths += 1
            if depth > float(min_depth):
                passing.append(fields[indv_i])
    if finite_depths == 0:
        passing = list(samples)
        with open(log_path, "a") as handle:
            handle.write(
                "No finite per-sample SNP depth values were available; "
                "skipping min_depth sample filter and retaining all VCF samples.\n"
            )
    with open(samps_path, "w") as handle:
        for sample in passing:
            handle.write(f"{sample}\n")


def _shell(args, log_path, append=False):
    mode = "ab" if append else "wb"
    with open(log_path, mode) as handle:
        subprocess.run(args, stdout=handle, stderr=handle, check=False)


def main():
    if len(sys.argv) != 10:
        print(__doc__)
        sys.exit(2)
    (vcf, prefix, min_depth, log_path,
     depth_out, miss_out, samps_out, summ_out, het_out) = sys.argv[1:]

    os.makedirs(os.path.dirname(depth_out), exist_ok=True)
    format_ids, samples = _read_vcf_header(vcf)

    _shell(["vcftools", "--gzvcf", vcf, "--FILTER-summary", "--out", prefix], log_path)
    if "DP" in format_ids:
        with open(log_path, "a") as handle:
            handle.write("FORMAT/DP found; computing per-sample SNP depth with vcftools --depth.\n")
        _shell(["vcftools", "--gzvcf", vcf, "--out", prefix, "--depth"], log_path, append=True)
        _normalize_depth_file(depth_out)
    elif "AD" in format_ids:
        with open(log_path, "a") as handle:
            handle.write("FORMAT/AD absent and FORMAT/AD found; computing per-sample SNP depth from summed AD values.\n")
        _write_ad_depth(vcf, depth_out, samples)
    else:
        with open(log_path, "a") as handle:
            handle.write("FORMAT/DP and FORMAT/AD absent; writing unavailable per-sample SNP depth values.\n")
        _write_empty_depth(depth_out, samples)

    _shell(["vcftools", "--gzvcf", vcf, "--out", prefix, "--het"], log_path, append=True)
    _shell(["vcftools", "--gzvcf", vcf, "--out", prefix, "--missing-indv"], log_path, append=True)
    _write_sample_filter(depth_out, samples, samps_out, min_depth, log_path)


if __name__ == "__main__":
    main()
