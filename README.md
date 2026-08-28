# oxo-flow-snparcher — Variant calling for non-model organisms: trimming, alignment and per-sample gVCFs

> ★ Verified · ⇄ Official port of [`harvardinformatics/snparcher`](https://github.com/harvardinformatics/snparcher) @ `v2.2` — same tools, same versions, same commands. Part of the [oxo-flow-community catalog](https://oxo-flow-community.github.io/).

[![CI](https://github.com/oxo-flow-community/oxo-flow-snparcher/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-snparcher/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

Point this workflow at a directory of paired-end FASTQ files plus a reference
genome, and it returns per-sample gVCFs, a cohort joint-genotyping VCF,
hard-filtered clean SNP/indel call sets, a callable-sites BED, and an
interactive QC dashboard. Each read pair is trimmed and filtered with fastp,
aligned to your reference with BWA-MEM, and called with GATK HaplotypeCaller
(or DeepVariant, opt-in) using low-coverage-appropriate defaults
(`-ploidy 2`, `--emit-ref-confidence GVCF`, `--min-pruning 1`); a QC stage
aggregates fastp and samtools metrics for every sample into a single TSV
report. The pipeline is designed for non-model organisms where no population
panel exists — it assumes diploid calls from a single sample at a time.

Input types match upstream's sample sheet: local paired FASTQs (default),
SRA accessions (downloaded with sra-tools), or external BAMs (staged into the
workflow). All upstream branches are ported except the structurally
impossible ones listed in the [Fidelity](#fidelity) table.

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
  1.24, gatk4 4.6.2.0, sambamba 1.0.1, sra-tools 3.2.1, mosdepth 0.3.3,
  clam, genmap, glnexus, vcftools 0.1.16, plink2/plink, admixture, R with
  tidyverse/plotly/flexdashboard), so conda or mamba must be installed at
  runtime; oxo-flow creates the environments on first run. The DeepVariant
  branch additionally requires a Docker-capable backend (runs
  `google/deepvariant:1.10.0`).

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
library_id, mark_duplicates`) is represented by the sample groups in
`main.oxoflow` (`[[sample_groups]]`). Each sample has one input unit
(empty `library_id` defaults to `sample_id`; multi-row sheets with multiple
libraries/units are not ported) and a metadata `input_type`:

| input_type | upstream `input` | oxo-flow declaration | branch |
|---|---|---|---|
| `fastq` (default) | local read pairs | files in `raw/<sample>_{1,2}.fastq.gz` | `fastp` → alignment |
| `srr` | SRA accession | group metadata `accession = "SRR..."` | `download_sra` (prefetch/fasterq-dump, ffq ENA fallback) → `fastp_srr` → alignment |
| `bam` | external BAM path | group metadata `bam_path = "/path/to/sample.bam"` | `stage_external_bam` (symlink) → indexing/calling |
| `gvcf` | external gVCF path | group metadata `gvcf_path = "/path/to/sample.g.vcf.gz"` | `normalize_external_gvcf_for_gatk` (recompress + index) → joint genotyping; gVCF samples skip trimming/alignment/calling entirely, exactly like upstream |

Example SRA, external-BAM and external-gVCF groups are documented in the
header of `main.oxoflow`. Each gvcf sample needs its own group (the
`gvcf_path` metadata is per group, like `bam_path`). The per-sample
`mark_duplicates` column is mapped to the global `mark_duplicates` config
key (upstream applies it per sample), with an optional per-sample override
via the `[workflow]` `metadata_file` table — see [Per-sample
`mark_duplicates` override](#per-sample-mark_duplicates-override).

### Per-sample `mark_duplicates` override

Upstream snpArcher reads `mark_duplicates` from each sample-sheet row
(default `true`); the port maps it to the global `config.mark_duplicates`
key (default `false`). To override the global key per sample, set
`[workflow]` `metadata_file` to a sample-keyed TSV/CSV/JSON table with a
`mark_duplicates` column:

```text
# metadata/samples.tsv
sample	mark_duplicates
sample1	true
sample2	false
```

```toml
[workflow]
metadata_file = "metadata/samples.tsv"
```

Precedence: a sample whose row defines the column uses that row's value
(`true` routes it through sambamba `markdup_library` →
`merge_dedup_libraries`; `false` through the no-markdup merge path). A
sample with no row — or a table without the column — falls back to the
global `config.mark_duplicates` exactly as before, so workflows that never
add the column behave byte-identically to earlier versions. Values must be
lowercase `true`/`false`; any other value is treated as "no override" and
falls back to the global key. The table is a lookup only: it does not
define the sample set (`[[sample_groups]]` does), and `input_type` routing
is unaffected. The `when` gates are the two-part form
`(config.mark_duplicates || {meta.mark_duplicates} == 'true') &&
{meta.mark_duplicates} != 'false'` (markdup) and
`(!config.mark_duplicates || {meta.mark_duplicates} == 'false') &&
{meta.mark_duplicates} != 'true'` (no-markdup) — each reduces exactly to
today's `config.mark_duplicates` truthiness when the column is absent.
This uses the engine's `metadata_file` + `{meta.*}` `when`-baking (engine
PR #234, post-v0.16.0); the CI builds the engine from source accordingly.

### Branch toggles

Every upstream branch is gated by a config key, all defaulting to the same
values as upstream's default config — with the port's documented deviation
that `mark_duplicates` and `joint_genotyping` default OFF so the default plan
runs the same 12 rules the port has always run (upstream defaults: both ON).

| config key | default | activates |
|---|---|---|
| `mark_duplicates` | `false` | sambamba `markdup_library` → `merge_dedup_libraries` (upstream default `true`); per-sample override via the optional `[workflow]` `metadata_file` `mark_duplicates` column |
| `variant_tool` | `"gatk"` | `"deepvariant"` switches per-sample calling to `deepvariant_call*` (docker); `"bcftools"` switches joint calling to the per-region `bcftools_call*` (see Fidelity) |
| `joint_genotyping_enabled` | `false` | GATK `create_db_mapfile` → `joint_genomics_db_import` → `joint_genotype_gvcfs`; for DeepVariant: `glnexus_joint` (upstream default `true`) |
| `generate_filtered_vcf` | `false` | `variant_filtration` GATK hard filters (upstream default `true`) |
| `callable_sites_enabled` | `false` | `mosdepth*`/clam coverage + genmap mappability → `callable_sites.bed` (upstream default `true`) |
| `modules_postprocess_enabled` | `false` | postprocess module (upstream default `false`; requires callable sites + joint genotyping) |
| `modules_qc_enabled` | `false` | qc module (upstream default `false`; requires joint genotyping — consumes the joint VCF) |
| `intervals_enabled` | `false` | interval scatter: `picard_intervals` → per-interval `gatk_haplotypecaller_interval*` → `concat_interval_gvcfs` (gVCF mode) and `create_db_intervals` → `gatk_genomics_db_import_interval` → `gatk_genotype_gvcfs_interval` → `concat_interval_vcfs` (joint mode). Runtime fan-out via engine `output_pattern`; requires engine ≥ 0.17 (upstream default `true` — see deviations) |
| `intervals_scatter_count` | `50` | `--scatter-count` of both SplitIntervals calls (upstream `num_gvcf_intervals`) |
| `intervals_db_scatter_factor` | `0.15` | DB shard count = `db_scatter_factor × samples × num_gvcf_intervals` (upstream `db_scatter_factor`) |
| `intervals_min_nmer` | `500` | picard `ScatterIntervalsByNs MIN_NMER` (upstream `min_nmer`) |
| `intervals_min_contig_length` | `0` | filter_picard_intervals `--min-contig-length` (upstream `min_contig_length`) |
| `intervals_db_max_intervals_per_shard` | `200` | split-db shard cap (upstream `db_max_intervals_per_shard`) |
| `intervals_db_max_contigs_per_shard` | `200` | split-db shard cap (upstream `db_max_contigs_per_shard`) |
| `bcftools_min_mapq` | `20` | bcftools mpileup `-q` (upstream `min_mapq`) |
| `bcftools_min_baseq` | `20` | bcftools mpileup `-Q` (upstream `min_baseq`) |
| `bcftools_max_depth` | `250` | bcftools mpileup `-d` (upstream `max_depth`) |

Parameter keys (`gatk_het_prior`, `deepvariant_model_type`,
`callable_sites_min_coverage`, `callable_sites_max_coverage`,
`callable_sites_fraction`, `callable_sites_merge_distance`,
`callable_sites_kmer`, `callable_sites_min_score`, `postprocess_contig_size`,
`postprocess_maf`, `postprocess_missingness`, `postprocess_exclude_scaffolds`,
`qc_min_depth`, `qc_max_sample_missingness`, `qc_exclude_scaffolds`,
`qc_clusters`, `qc_google_api_key`, `sample_metadata`) mirror upstream's
config values; pass overrides as `oxo-flow run main.oxoflow key=value`.
`sample_metadata` is upstream's optional qc-module metadata CSV
(sample_id/long/lat columns) that feeds the dashboard's terrain-map panel —
set it alongside `modules_qc_enabled=true` and `qc_google_api_key`.

## Live verification

- `normalize_external_gvcf_for_gatk` + `joint_genomics_db_import` +
  `joint_genotype_gvcfs` + `generate_coords_file`: live-verified on
  bioinfo-wsx 2026-08-28 (external gVCF cohort, 0 failed; the run also
  fixed two gvcf-cohort bugs — `parse_bam_stats` and
  `combine_qc_metrics.py` now skip gvcf samples, PR #13).
- bcftools caller branch (`bcftools_regions` → `bcftools_call` →
  `bcftools_concat_regions`): live-verified 2026-08-28 on macOS (engine
  0.16.0 + output_pattern from engine-235, conda). The runtime `.fai`
  scan instantiated one deferred `bcftools_call` per contig (1-contig
  fixture → exactly one region) and produced real variants in
  `results/vcfs/regions/L000000.vcf.gz` (GT:PL:AD for all samples,
  AC/AN=6); 30 succeeded / 0 failed / 34 skipped.
- interval-scatter branch (`picard_intervals` → `create_gvcf_intervals` /
  `create_db_intervals` → per-interval HaplotypeCallers / per-shard
  GenomicsDBImport + GenotypeGVCFs → `concat_interval_gvcfs` /
  `concat_interval_vcfs`): live-verified 2026-08-28 on macOS (same
  engine). The SplitIntervals scans instantiated per-interval
  HaplotypeCallers per sample (1-contig fixture → 1 interval) and 22 db
  shards (DB_SCATTER = 0.15 × 3 samples × 50), each shard run through
  GenomicsDBImport + GenotypeGVCFs and concatenated to
  `results/vcfs/raw.vcf.gz` (30 real variants, all 3 samples); 85
  succeeded / 0 failed / 37 skipped, 117 outputs verified.
- The qc module's plink steps need **≥ 2 samples** (a 1-sample cohort
  prunes to an empty VCF and `qc_plink` fails with "No samples in .vcf
  file" — upstream behaves the same on a single-sample cohort).

## Source

Ported from **[harvardinformatics/snparcher](https://github.com/harvardinformatics/snparcher)**,
version `v2.2` (commit `e0e7a9478d4e042fce217db4e6077dafdaf57245`, MIT).
Created 2026-08-15; this workflow may lag behind upstream releases. Upstream
license and attribution are recorded in [NOTICE.md](NOTICE.md).

## Fidelity

74 rules ported from upstream v2.2 (up from 60), covering every branch that
can be expressed in oxo-flow's DAG. Commands are ported verbatim
(same flags, same output paths); upstream's snakemake `{{...}}` shell escaping
is unwrapped to literal braces, and snakemake `{params.*}`/`{resources.*}`
references are resolved to their upstream values or to `{config.*}` keys.
The two runtime-fan-out branches (interval scatter, bcftools caller) are
ported with the engine's `output_pattern` primitive (issue #227 item 5,
merged in oxo-flow 0.17): a producer rule whose outputs are enumerated by a
filesystem scan after it completes, instantiating one downstream consumer per
discovered value. These branches require an engine with `output_pattern`
support; older engines ignore the key and the when-gates keep the default
plan unchanged.

| Upstream rule | oxo-flow rule | Tool (version) | Notes |
|---|---|---|---|
| `prepare_reference` (local branch) | `prepare_reference` | samtools 1.24 (bgzip) | identical command; url/accession branches not ported (config `reference_source` is a local path) |
| `index_reference` | `index_reference` | samtools 1.24, bwa 0.7.19 | identical command (faidx + dict + bwa index) |
| `fastp` | `fastp` / `fastp_srr` | fastp 1.3.6 | identical flags; `{sample}/{sample}/u1` fan-out for single-row sheets with empty `library_id`; SRA-downloaded reads via `fastp_srr` |
| `download_sra` | `download_sra` | sra-tools 3.2.1, ffq, curl, pigz | identical prefetch→fasterq-dump→pigz flow with ffq/ENA fallback; fasterq-dump `--tmpdir` dropped (oxo-flow has no per-rule tmpdir) |
| `bwa_mem` | `bwa_mem` | bwa 0.7.19, samtools 1.24 | identical command incl. read group `ID:{sample}.u1 SM:{sample} LB:{sample} PL:ILLUMINA`; raw BAM is temp like upstream |
| `merge_library_bams` | `merge_library_bams` | samtools 1.24 | per-library merge, single input unit in the default path |
| `merge_library_level_bams` | `merge_library_level_bams` | samtools 1.24 | no-markdup path (`results/bams/merged/{sample}.bam`) |
| `markdup_library` / `merge_dedup_libraries` | `markdup_library` / `merge_dedup_libraries` | sambamba 1.0.1, samtools 1.24 | identical commands; gated on `mark_duplicates` (default `false`; upstream default `true` — see deviations below; per-sample override via the `metadata_file` `mark_duplicates` column) |
| `index_bam_csi` | `index_bam_csi` / `index_bam_csi_markdup` / `index_bam_csi_external` | samtools 1.24 | identical (`samtools index -c`), one per BAM-producing branch |
| `stage_external_bam` | `stage_external_bam` | — | external BAM inputs symlinked into `results/bams/input/` then run through the standard callers |
| `normalize_external_gvcf_for_gatk` / `archive_gatk_gvcf` (gvcf input type) | `normalize_external_gvcf_for_gatk` | bcftools 1.23 | external gVCF inputs recompressed + tabix-indexed to `results/gvcfs/{sample}.g.vcf.gz` (upstream long-contig mode's archive command) and fed straight into joint genotyping; gVCF samples skip calling. Upstream short mode feeds the raw external path to the mapfile; the port normalizes so the uniform `results/gvcfs/{sample}.g.vcf.gz` pattern holds. Upstream refuses gvcf inputs with non-GATK callers; the port accepts them in the GLnexus path (normalized gVCFs are valid GLnexus input) — see deviations |
| `gatk_haplotypecaller` (standard mode) | `gatk_haplotypecaller` / `_markdup` / `_external` | gatk4 4.6.2.0 | identical flags incl. `-ploidy 2 --emit-ref-confidence GVCF --min-pruning 1 --min-dangling-branch-length 1` (low-coverage defaults); `-Xmx7000m` = upstream default profile `mem_mb_reduced`; threads 1 as upstream |
| `picard_intervals` (interval mode) | `picard_intervals` | picard 3.5.0 | identical `ScatterIntervalsByNs` call (`MAX_TO_MERGE=500`, `OUTPUT_TYPE=ACGT` = upstream `ScatterIntervalsByNs` block) |
| `create_gvcf_intervals` (interval mode) | `create_gvcf_intervals` | gatk4 4.6.2.0 | upstream `intervals.smk:59` checkpoint: one `SplitIntervals --scatter-count N --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION` per sample, then the per-interval file list is enumerated by the engine's `output_pattern` scan (see deviations) |
| `gatk_haplotypecaller` (interval mode) | `gatk_haplotypecaller_interval` / `_markdup` / `_external` (short mode) + `_long` twins | gatk4 4.6.2.0 | deferred consumers, instantiated once per discovered interval; same flags as standard mode plus `-L {interval}` from the scattered list; long-contig mode writes plain `.vcf` (no `.gz`, no index) since GATK silently omits the index for `-O *.vcf.gz` on long contigs |
| `concat_interval_gvcfs` | `concat_interval_gvcfs` | bcftools 1.23 | `bcftools concat -D -a` over the per-interval gVCFs + `bcftools sort` + `index -t`, one per sample |
| `create_db_intervals` (interval mode) | `create_db_intervals` | gatk4 4.6.2.0, python | upstream `intervals.smk:89` checkpoint: `DB_SCATTER = db_scatter_factor × samples × num_gvcf_intervals` (computed from `{config.samples_list}`), `SplitIntervals --subdivision-mode INTERVAL_SUBDIVISION`, then `scripts/interval_list_tools.py split-db` re-shards to the per-shard interval/contig caps |
| `gatk_genomics_db_import` (interval mode) | `gatk_genomics_db_import_interval` | gatk4 4.6.2.0 | deferred consumer, one import per db shard (`-L {db_interval}`, mapfile + `--merge-input-intervals` like the cohort rule) |
| `gatk_genotype_gvcfs` (interval mode) | `gatk_genotype_gvcfs_interval` | gatk4 4.6.2.0 | deferred consumer, `gendb://` per-shard GenotypeGVCFs → `results/vcfs/intervals/L{db_interval}.vcf.gz` |
| `concat_interval_vcfs` | `concat_interval_vcfs` / `concat_interval_vcfs_long` | bcftools 1.23 | `bcftools concat -D -a` + sort + index over the per-shard VCFs → `results/vcfs/raw.vcf.gz` (+ `.tbi` short mode, `.csi` long mode) |
| `bcftools_regions` (bcftools caller) | `bcftools_regions` | bcftools 1.23 | upstream `bcftools.smk:13` checkpoint: enumerates reference contigs from the runtime `.fai` into per-contig marker files + `regions.tsv`; the file list is runtime-discovered via `output_pattern` |
| `bcftools_call` (bcftools caller) | `bcftools_call` | bcftools 1.23 | deferred consumer, one mpileup+call per region (`-q/-Q/-d` from config, `-r "$CONTIG"`, `--ploidy`, `-v`); the cohort BAM list is resolved per sample with upstream's markdup → merged → external priority (see deviations) |
| `bcftools_concatenate_vcfs` (bcftools caller) | `bcftools_concat_regions` | bcftools 1.23 | `bcftools concat -D -a` + sort + index over the per-region VCFs → `results/vcfs/raw.vcf.gz` |
| `deepvariant_call` | `deepvariant_call` / `_markdup` / `_external` | google/deepvariant:1.10.0 (docker) | identical `/opt/deepvariant/bin/run_deepvariant` invocation; gated on `variant_tool = "deepvariant"` |
| `create_db_mapfile` | `create_db_mapfile` | python (script) | identical logic, ported as `scripts/write_joint_gvcf_mapfile.py` |
| `joint_genomics_db_import` | `joint_genomics_db_import` | gatk4 | identical GenomicsDBImport flow incl. `TILEDB_DISABLE_FILE_LOCKING` and `--merge-input-intervals` from `scripts/interval_list_tools.py` (merge threshold 50 = upstream `GENOMICSDB_MERGE_CONTIG_THRESHOLD`) |
| `joint_genotype_gvcfs` | `joint_genotype_gvcfs` | gatk4 | identical (tar-extract → `gendb://` GenotypeGVCFs → `results/vcfs/raw.vcf.gz`); temp raw VCF like upstream |
| `glnexus_joint` | `glnexus_joint` | glnexus, bcftools 1.23 | identical DeepVariant-config GLnexus join; `mem_gbytes` = `mem_mb_reduced/1024` rounded to 8, computed from the default profile (see deviations) |
| `variant_filtration` | `variant_filtration` | gatk4, bcftools | identical RPRS/FS_SOR/MQ/QUAL hard filters, `--invalidate-previous-filters true`, then `bcftools index -f -t` |
| long-contig (CSI) indexing — `_resolve_long_contig_mode` (common.smk) + `compress_interval_raw_vcf` (interval mode) + postprocess module `regions_to_index` | `resolve_long_contig_mode` + `*_long` twins (`gatk_haplotypecaller_interval_long`/`_markdup_long`/`_external_long`, `concat_interval_gvcfs_long`, `create_db_mapfile_long`, `gatk_genomics_db_import_interval_long`, `gatk_genotype_gvcfs_interval_long`, `concat_interval_vcfs_long`, `variant_filtration_long`, `postprocess_basic_filter_long`, `postprocess_strict_filter_long`, `postprocess_subset_indels_long`, `postprocess_subset_snps_long`, `postprocess_drop_indel_snps_long`) | bcftools 1.24, gatk4 4.6.2.0 | config `long_contig_mode` ("auto"/true/false) mirrors upstream `TBI_MAX_CONTIG_LENGTH = 2**29 - 1` auto-detection from the `.fai`; long mode emits `.csi` (`bcftools index -c`) and keeps GATK consumers on plain `.vcf` + `.idx`; the postprocess module follows the same twin pattern (`.csi` produced and consumed end-to-end); short mode is byte-identical to a run without the key (see deviations 3, 10, 13, 14) |
| `collect_fastp_stats` | `collect_fastp_stats` | python (script) | identical logic, ported as `scripts/collect_fastp_stats.py` |
| `bam_stats` | `bam_stats` / `_markdup` / `_external` | samtools 1.24 | identical (coverage + flagstat -O tsv); outputs temp like upstream |
| `parse_bam_stats` | `parse_bam_stats` | python (script) | identical logic, ported as `scripts/parse_bam_stats.py` |
| `combine_qc_metrics` | `combine_qc_metrics` | python (script) | identical report format; gather via `expand_inputs` |
| `mosdepth` | `mosdepth` / `_markdup` / `_external` | mosdepth 0.3.3 | identical (`--d4 -t {threads}`), per BAM branch |
| `clam_collect` | `clam_collect` | clam | identical (`clam collect -o depths.zarr`) |
| `callable_coverage_thresholds` | `callable_coverage_thresholds` | python (script) | identical logic, ported verbatim as `scripts/callable_coverage_thresholds.py` |
| `clam_loci` | `clam_loci` | clam | identical incl. per-sample mode and `-m/-M` from the thresholds TSV |
| `coverage_bed` | `coverage_bed` | python, bedtools | identical logic, ported verbatim as `scripts/callable_zarr_to_bed.py` |
| `genmap_index` | `genmap_index` | genmap | identical index-mode switching on decompressed FASTA size (skew ≥ 5 GB, sampled ≥ 2 GB) |
| `genmap_mappability` | `genmap_mappability` | genmap | identical (`-K 150 -E 2 -bg -T`) |
| `mappability_bed` | `mappability_bed` | awk, bedtools | identical score filter + `-d 100` merge |
| `callable_sites_bed` | `callable_sites_bed` | bedtools | identical sort/merge of coverage + mappability BEDs |
| postprocess module (`filter_individuals`, `basic_filter`, `update_bed`, `strict_filter`, `subset_snps`, `subset_indels`, `drop_indel_SNPs`) | `postprocess_filter_individuals` … `postprocess_drop_indel_snps` | bcftools 1.23, awk, bedtools, tabix | identical commands; `scripts/write_include_samples.py` for the sample list; AF upper bound computed as `1 - maf` (upstream `1-{params.maf}`) |
| qc module (`contig_map`, `vcftools_individuals`, `subsample_snps`, `prepare_plink_inputs`, `copy_qc_report`, `plink`, `setup_admixture`, `admixture`, `generate_coords_file`, `qc_dashboard`) | `qc_contig_map` … `qc_dashboard` | vcftools 0.1.16, bcftools 1.23, plink2/plink, admixture, R | identical commands; logic ported to `scripts/contig_map.py`, `vcftools_individuals.py`, `prepare_plink_inputs.py`, `contigs4admixture.py`, `generate_coords.py`, `qc_dashboard_render.R`. `generate_coords_file` reads the optional `sample_metadata` CSV (lat/long) into `results/qc/coords.txt`, consumed by the dashboard's terrain-map panel when `qc_google_api_key` is set; without metadata it writes the empty placeholder (upstream's own else branch) so the panel just prints its placeholder text |
| `setup` / `download_reads` / `map_samples` / `call_variants` / `qc_report` / `callable_sites` / `gvcfs` (Snakefile aggregation targets) | n/a | — | Snakemake target rules, no commands of their own |

### Remaining exclusions

| Item | Why excluded | Evidence |
|---|---|---|
| parabricks (all `parabricks_*` rules) | requires `--nv` GPU passthrough (upstream `parabricks.smk` runs `--nv` images with `nvidia-docker`); the oxo-flow docker backend has no `--nv` support and no GPU device declaration; additionally NVIDIA EULA/license enforcement cannot be guaranteed in CI | `variant_calling/parabricks.smk` (every rule is `--nv`) |
| sentieon (all `sentieon_*` rules) | proprietary tool gated on a `SENTIEON_LICENSE` server and a pre-installed license; cannot be distributed or verified in a community port | `config/config.yaml` `sentieon` section; `workflow/rules/sentieon.smk` |
| `denovo` and `structural_variants` pipeline sections | do not exist as rules in upstream v2.2 | grep of `workflow/` at e0e7a94 finds neither rule set |
| multi-library / multi-unit rows (library_id, input_unit) | the sample-group model has one unit per sample; consumers of `results/bams/raw/{sample}/{sample}/u1.bam` are hard-coded to the `u1` unit | sample sheet semantics in upstream `README` |

### Documented deviations from upstream

1. **Default config differs from upstream**: upstream defaults `mark_duplicates: true`, `joint_genotyping: true`, `generate_filtered_vcf: true`, `callable_sites.enabled: true`; the port defaults all of these to `false` so its default plan is byte-identical to the previous 12-rule port. Flip the keys to get upstream's full pipeline.
2. **postprocess/qc modules consume `results/vcfs/raw.vcf.gz`**, not upstream's `FINAL_VCF` (which is the hard-filtered VCF when `generate_filtered_vcf: true` and GATK is used). The modules were run on the raw joint VCF. The difference only matters when combining GATK + `generate_filtered_vcf` + a module, and reproduces upstream behavior for the DeepVariant path.
3. **Postprocess module long-contig CSI mode**: upstream conditionally uses CSI indexes when contigs exceed 512 Mb (`regions_to_index`); the port's postprocess rules follow the calling side's `*_long` twin pattern — `postprocess_basic_filter_long`/`postprocess_strict_filter_long`/`postprocess_subset_indels_long`/`postprocess_subset_snps_long`/`postprocess_drop_indel_snps_long` emit `.csi` (`bcftools index -c`; `tabix -C` for the SNP-positions file) and consume the `.csi`-indexed upstream VCFs of the chain when `long_contig_mode` selects long mode. Short mode is unchanged (`.tbi` everywhere, byte-identical plans).
4. **`glnexus_joint` memory**: upstream computes `mem_gbytes` from the default profile's `mem_mb_reduced`; the port inlines the resulting value 8 (with the same `if < 1 then 1` clamp).
5. **QC-metrics and callable-sites branches with non-fastq input types**: `combine_qc_metrics` and the callable-sites `expand_inputs` reference `results/fastp/{sample}/{sample}/u1.json` / `results/callable_sites/depths/{sample}.*` for every sample, which only exist for fastq/srr (fastp stats) or BAM-bearing samples (depths). For bam/gvcf cohorts the expands fail at plan time. Use fastq/srr groups when combining QC metrics or computing callable sites.
6. **`fasterq-dump --tmpdir` dropped** (no per-rule tmpdir in oxo-flow); SRA downloads use the current directory.
7. **gvcf inputs are accepted on any caller**: upstream hard-fails gvcf inputs with non-GATK callers (bcftools/deepvariant/parabricks); the port normalizes them regardless, so the DeepVariant GLnexus path also accepts them. Normalized gVCFs are valid GLnexus input, so this is a relaxation, not a behavior change.
8. **`coords.txt` is always produced in qc mode**: upstream only creates it when the metadata CSV actually has lat/long rows; the port writes the same file empty otherwise (upstream's own placeholder branch), so the dashboard's map panel shows its placeholder text instead of being absent.
9. **Interval scatter runs one `SplitIntervals` per sample** (gVCF mode): upstream's `create_gvcf_intervals` checkpoint runs once for the whole cohort and downstream per-interval rules address the shared `results/intervals/{sample}/...` per-sample subdirectory; the engine's `output_pattern` gives each producer instance its own scan domain, so the port runs the split per sample (identical inputs → identical interval lists) and the per-interval HaplotypeCallers read `results/intervals/gvcf/{sample}/{interval}-scattered.interval_list`. Functionally equivalent, N × the split work (N = samples); the per-sample scans are independent, so a scatter of 50 × 20 samples costs 20 short GATK calls instead of 1.
10. **No standalone `compress_interval_raw_vcf` step**: upstream re-compresses the per-shard raw VCFs before concatenation; the port's `concat_interval_vcfs` consumes the per-shard VCFs directly (`bcftools concat -D -a` is format-agnostic) and `bcftools sort` normalizes coordinate order. The long-contig (CSI) side of that step IS ported (`concat_interval_vcfs_long` emits `raw.vcf.gz` + `.csi`; `archive_gatk_gvcf` is folded into `concat_interval_gvcfs_long`, which writes the GATK-consumable plain `.vcf` + `.idx` that upstream stages in `results/gvcfs/work/` and never emits a separate durable gz+csi archive — nothing consumes one in the port's interval chain).
11. **`bcftools_call` resolves the per-sample final BAM in the shell** (upstream's `get_final_bam` markdup → merged → external priority) from `{config.samples_list}`, instead of a static input per branch; the declared `results/bams/{merged,markdup,input}/*.bam` glob inputs still order the calls behind every BAM-producing branch via DAG edges. A sample with no final BAM (e.g. gvcf-only cohorts) contributes nothing and the call proceeds with the remaining samples; an empty BAM list exits 0 with a log note (upstream's shell has the same degenerate-cohort behavior).
12. **The interval and bcftools branches need an engine with `output_pattern`** (oxo-flow ≥ 0.17). On older engines the key is ignored and the fresh wildcards (`{interval}`, `{db_interval}`, `{region}`) stay unbound, so `validate`/`lint`/`dry-run` still pass (when-gates keep both branches off by default) but those branches must not be enabled there.
13. **Long-contig mode gates on a side-effect flag file**: `resolve_long_contig_mode` writes a fixed-path flag (`results/reference/.long_contig.flag`) and every `*_long` rule is gated on it (`depends_on = ["resolve_long_contig_mode"]`); the short rules carry the negated gate. Engine gotcha: a `when` containing `wildcard.` is evaluated at PLAN time for DAG morphing, when the flag does not exist yet — a plan-level long rule would be dropped from the plan entirely. Plan-level long rules therefore avoid `wildcard.` in `when`: `concat_interval_gvcfs_long` gates on `{meta.input_type} != 'gvcf'` (baked per instance at expansion time, no `wildcard.` reference), so its `file_exists(...)` gate is evaluated at execution time. Rules instantiated at run time via `output_pattern` (per-interval HC long, db import long, genotype long) are unaffected. The `when` evaluator does not expand `{config.x}` inside `file_exists()`, hence the fixed flag path. Changing the reference or `long_contig_mode` between runs requires a fresh run directory (standard checkpoint semantics).
14. **GATK cannot read `.csi`-indexed VCFs** (verified with GATK 4.6.2.0: VariantFiltration and GenomicsDBImport reject `.csi`, accepting `.tbi`/`.idx`; GATK also silently writes no index for `-O x.vcf.gz` on long contigs). Upstream's long-mode `variant_filtration` feeds `raw.vcf.gz` + `.csi` to GATK and fails; the port's `variant_filtration_long` converts to plain `.vcf` in-shell (`bcftools view -O v`), filters, and re-indexes with `bcftools index -f -c`. GenomicsDBImport reads the plain work gVCFs via GATK `.idx` (created by `gatk IndexFeatureFile` in `concat_interval_gvcfs_long`); GATK reads plain `.vcf` without any index for streaming tools, but db import strictly requires `.idx`. Memory: the port mirrors upstream's default GATK heap (`-Xmx7000m` = upstream profile `mem_mb_reduced: attempt * 7000`). A fully homozygous-reference interval over a contig at the CSI limit can exhaust this heap at gVCF finalization (GATK materializes the whole hom-ref block's depth array for `getMedianDP`; a 600 Mb hom-ref block peaks around 7.2 GB and dies with `OutOfMemoryError` at the default heap — verified live, upstream-inherited since upstream assigns whole contigs to intervals too, `BALANCING_WITHOUT_INTERVAL_SUBDIVISION`). Raise `-Xmx` in the two interval HaplotypeCaller rules (short and long) for such references. The postprocess module is bcftools-only, so its long mode needs no GATK workaround: the `*_long` twins produce and consume `.csi` directly. Two long-chain fixes on top of PR #17: `concat_interval_gvcfs_long` and `concat_interval_vcfs_long` stage bgzip-compressed, indexed copies of the plain per-interval gVCFs / per-shard VCFs in `results/{gvcfs,vcfs}/work/` before concatenating — `bcftools concat` cannot read the plain `.g.vcf`/`.vcf` files the long HaplotypeCaller and GenotypeGVCFs emit (GATK `.idx` convention). One long-chain fix on top of PR #17: `concat_interval_gvcfs_long` stages bgzip-compressed, indexed copies of the plain per-interval gVCFs in `results/gvcfs/work/` before concatenating — `bcftools concat` cannot read the plain `.g.vcf` files the long HaplotypeCaller emits.

Version pinning: upstream envs declare only `>=` ranges with no lockfile;
exact pins (fastp 1.3.6, samtools 1.24, bwa 0.7.19, gatk4 4.6.2.0, bcftools
1.23, sra-tools 3.2.1, mosdepth 0.3.3, vcftools 0.1.16, plink2) were resolved
from bioconda/conda-forge at port time (2026-08-15). Upstream default-profile
thread overrides (fastp 6, bwa_mem 16) are runtime knobs; the port keeps the
rules' own declarations (4 and 8).

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
