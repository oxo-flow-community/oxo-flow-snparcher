# oxo-flow-snparcher — Variant calling for non-model organisms: trimming, alignment and per-sample gVCFs

> ★ Verified · ⇄ Official port of [`harvardinformatics/snparcher`](https://github.com/harvardinformatics/snparcher) @ `v2.2` — same tools, same versions, same commands. Part of the [oxo-flow-community catalog](https://oxo-flow-community.github.io/).

[![CI](https://github.com/oxo-flow-community/oxo-flow-snparcher/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-snparcher/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

Point this workflow at a directory of paired-end FASTQ files plus a reference
genome, and it returns per-sample gVCFs and a cohort QC report. Each read pair
is trimmed and filtered with fastp, aligned to your reference with BWA-MEM,
and called with GATK HaplotypeCaller using low-coverage-appropriate defaults
(`-ploidy 2`, `--emit-ref-confidence GVCF`, `--min-pruning 1`); a QC stage
aggregates fastp and samtools metrics for every sample into a single TSV
report. The pipeline is designed for non-model organisms where no population
panel exists — it assumes diploid calls from a single sample at a time.

## Installation

### 1. Install oxo-flow

Requires oxo-flow >= 0.12.0. Prebuilt release binary (recommended):

```bash
curl -fL -o oxo-flow.tar.gz \
  https://github.com/Traitome/oxo-flow/releases/latest/download/oxo-flow-latest-x86_64-unknown-linux-gnu.tar.gz
tar xzf oxo-flow.tar.gz
sudo mv oxo-flow /usr/local/bin/
```

Alternatively via conda: `conda install -c bioconda oxo-flow-cli` (note: the
bioconda package may lag behind releases; other platform binaries are on the
[releases page](https://github.com/Traitome/oxo-flow/releases)).

### 2. Get this workflow

```bash
git clone https://github.com/oxo-flow-community/oxo-flow-snparcher.git
cd oxo-flow-snparcher
```

### 3. Requirements

- **Reference data**: the only external reference input is your genome — a
  FASTA file (plain or gzip-compressed), passed at run time as
  `reference_source=/path/to/genome.fa.gz` (the committed default points at
  the tiny `test/fixtures/ref/genome.fa` fixture, so set it before real
  runs). The workflow bgzip-compresses and indexes it itself
  (`prepare_reference` + `index_reference`), so no pre-built indices are
  required. Reads go in `raw/<sample>_1.fastq.gz` / `raw/<sample>_2.fastq.gz`
  for every sample listed in the `[[sample_groups]]` of `main.oxoflow`.
- **Compute**: up to 8 CPUs and 8 GB per rule — `bwa_mem` uses 8 threads,
  `fastp` 4, and `gatk_haplotypecaller` a 7 GB Java heap (1 thread).
- **Tools**: delivered as conda environments with pinned versions
  (`envs/*.yaml`, conda-forge + bioconda: fastp 1.3.6, bwa 0.7.19, samtools
  1.24, gatk4 4.6.2.0), so conda or mamba must be installed at runtime;
  oxo-flow creates the environments on first run.

## Usage

```bash
# 1. install oxo-flow (see Requirements)
# 2. prepare data: raw/<sample>_1.fastq.gz / raw/<sample>_2.fastq.gz
#    (tiny examples in test/fixtures/raw/)
# 3. set your reference genome
oxo-flow run main.oxoflow reference_source=/path/to/genome.fa.gz
# 4. preview the plan
oxo-flow dry-run main.oxoflow
# 5. run a subset
oxo-flow run main.oxoflow --samples first:2
```

The sample sheet of upstream snpArcher (`sample_id, input_type, input,
library_id, mark_duplicates`) is represented by the sample list in
`main.oxoflow` (`[[sample_groups]]`) — one row per sample with local FASTQ
inputs, empty `library_id` (defaults to `sample_id`), and a single input
unit. Multi-row sheets (multiple libraries/units) and `srr`/`bam`/`gvcf`
input types are not ported.

## Source

Ported from **[harvardinformatics/snparcher](https://github.com/harvardinformatics/snparcher)**,
version `v2.2` (commit `e0e7a9478d4e042fce217db4e6077dafdaf57245`, MIT).
Created 2026-08-15; this workflow may lag behind upstream releases. Upstream
license and attribution are recorded in [NOTICE.md](NOTICE.md).

## Fidelity

Port scope: the default-parameters main execution path (FASTQ inputs,
`variant_calling.tool = gatk`), with the committee-approved exclusions
`markdup`, `joint_genotyping`, `denovo`, and `structural_variants` (the
latter two do not exist as steps in upstream v2.2). Intermediate machinery
that only serves excluded branches (interval scatter, GenomicsDB,
callable-sites, SRA downloads, optional modules, non-GATK callers) is listed
as "not ported" below.

| Upstream rule | oxo-flow rule | Tool (version) | Notes |
|---|---|---|---|
| `prepare_reference` (local branch) | `prepare_reference` | samtools 1.24 (bgzip) | identical command; url/accession branches not ported (config `reference_source` is a local path) |
| `index_reference` | `index_reference` | samtools 1.24, bwa 0.7.19 | identical command (faidx + dict + bwa index) |
| `fastp` | `fastp` | fastp 1.3.6 | identical flags; per (sample, library, input_unit) fan-out modeled as `{sample}/{sample}/u1` — the upstream default for single-row-per-sample sheets with empty `library_id` |
| `bwa_mem` | `bwa_mem` | bwa 0.7.19, samtools 1.24 | identical command incl. read group `ID:{sample}.u1 SM:{sample} LB:{sample} PL:ILLUMINA`; raw BAM is temp like upstream |
| `merge_library_bams` | `merge_library_bams` | samtools 1.24 | per-library merge, single input unit in the default path |
| `merge_library_level_bams` | `merge_library_level_bams` | samtools 1.24 | no-markdup path (`results/bams/merged/{sample}.bam`); used because markdup is excluded |
| `markdup_library` / `merge_dedup_libraries` | not ported | sambamba 1.0.1 | committee exclusion `markdup`; port default `mark_duplicates = false` (upstream default true) |
| `index_bam_csi` | `index_bam_csi` | samtools 1.24 | identical (`samtools index -c`) |
| `gatk_haplotypecaller` (standard mode) | `gatk_haplotypecaller` | gatk4 4.6.2.0 | identical flags incl. `-ploidy 2 --emit-ref-confidence GVCF --min-pruning 1 --min-dangling-branch-length 1` (low-coverage defaults); `-Xmx7000m` = upstream default profile `mem_mb_reduced`; threads 1 as upstream |
| `collect_fastp_stats` | `collect_fastp_stats` | python (script) | identical logic, ported as `scripts/collect_fastp_stats.py` |
| `bam_stats` | `bam_stats` | samtools 1.24 | identical (coverage + flagstat -O tsv); outputs temp like upstream |
| `parse_bam_stats` | `parse_bam_stats` | python (script) | identical logic, ported as `scripts/parse_bam_stats.py` |
| `combine_qc_metrics` | `combine_qc_metrics` | python (script) | identical report format; gather via `expand_inputs` |
| `download_sra` | not ported | sra-tools | SRA (`srr`) inputs out of scope; port default path uses local FASTQ inputs |
| `stage_external_bam` | not ported | — | BAM input type out of scope |
| interval machinery (`picard_intervals`, `create_gvcf_intervals`, `create_db_intervals`, `gatk_haplotypecaller_interval`, `concat_interval_gvcfs*`, `concat_interval_vcfs*`, `compress_interval_raw_vcf`, `normalize_external_gvcf_for_gatk`, `archive_gatk_gvcf`) | not ported | — | upstream default `intervals.enabled: true`; excluded per committee scope (interval scatter serves joint genotyping). Port default `intervals_enabled = false` |
| joint genotyping (`joint_genomics_db_import`, `joint_genotype_gvcfs`, `create_db_mapfile`, `gatk_genomics_db_import`, `gatk_genotype_gvcfs`) | not ported | — | committee exclusion `joint_genotyping`; per-sample gVCF is the port's final call-set output |
| `variant_filtration` (hard filters) | not ported | gatk4, bcftools | downstream of the excluded joint-genotyping raw VCF; cannot run faithfully without it |
| callable sites (`mosdepth`, `clam_collect`, `callable_coverage_thresholds`, `clam_loci`, `coverage_bed`, `genmap_index`, `genmap_mappability`, `mappability_bed`, `callable_sites_bed`) | not ported | mosdepth, clam, genmap, bedtools | coverage/mappability BED branch out of scope (not in the ~9-rule committee scope) |
| `bcftools_call` / `deepvariant_call` / `glnexus_joint` / parabricks / sentieon rules | not ported | — | non-GATK callers, selected by config only |
| postprocess module rules (`basic_filter`, `strict_filter`, `drop_indel_SNPs`, `subset_snps`, `subset_indels`, `contig_map`, `update_bed`, …) | not ported | — | `modules.postprocess.enabled: false` by default upstream |
| qc module rules (`plink`, `admixture`, `subsample_snps`, `filter_individuals`, `vcftools_individuals`, `prepare_plink_inputs`, `setup_admixture`, `generate_coords_file`, `copy_qc_report`, `qc_dashboard`, `denovo`-style dashboard inputs) | not ported | — | `modules.qc.enabled: false` by default upstream |
| `setup` / `download_reads` / `map_samples` / `call_variants` / `qc_report` / `callable_sites` / `gvcfs` (Snakefile aggregation targets) | n/a | — | Snakemake target rules, no commands of their own |

Version pinning: upstream envs declare only `>=` ranges with no lockfile;
exact pins (fastp 1.3.6, samtools 1.24, bwa 0.7.19, gatk4 4.6.2.0, picard
3.5.0, bcftools 1.24) were resolved from bioconda/conda-forge at port time
(2026-08-15). Upstream default-profile thread overrides (fastp 6, bwa_mem
16) are runtime knobs; the port keeps the rules' own declarations (4 and 8).

## Test

Run the acceptance suite (validate + lint + dry-run) against the committed
fixture data:

```bash
bash test/run.sh
```

## License

Apache-2.0. Copyright (c) 2026 oxo-flow-community. Upstream attribution in
[NOTICE.md](NOTICE.md). The upstream snpArcher project is MIT-licensed; its
license text is included verbatim at [LICENSE.upstream](LICENSE.upstream).
