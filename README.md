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
| `gvcf` | external gVCF path | not ported | see Fidelity — external paths cannot feed group expansion in a static DAG |

Example SRA and external-BAM groups are documented in the header of
`main.oxoflow`. The per-sample `mark_duplicates` column is mapped to the
global `mark_duplicates` config key (upstream applies it per sample); a
single group-level override is not ported.

### Branch toggles

Every upstream branch is gated by a config key, all defaulting to the same
values as upstream's default config — with the port's documented deviation
that `mark_duplicates` and `joint_genotyping` default OFF so the default plan
runs the same 12 rules the port has always run (upstream defaults: both ON).

| config key | default | activates |
|---|---|---|
| `mark_duplicates` | `false` | sambamba `markdup_library` → `merge_dedup_libraries` (upstream default `true`) |
| `variant_tool` | `"gatk"` | `"deepvariant"` switches per-sample calling to `deepvariant_call*` (docker) |
| `joint_genotyping_enabled` | `false` | GATK `create_db_mapfile` → `joint_genomics_db_import` → `joint_genotype_gvcfs`; for DeepVariant: `glnexus_joint` (upstream default `true`) |
| `generate_filtered_vcf` | `false` | `variant_filtration` GATK hard filters (upstream default `true`) |
| `callable_sites_enabled` | `false` | `mosdepth*`/clam coverage + genmap mappability → `callable_sites.bed` (upstream default `true`) |
| `modules_postprocess_enabled` | `false` | postprocess module (upstream default `false`; requires callable sites) |
| `modules_qc_enabled` | `false` | qc module (upstream default `false`) |
| `intervals_enabled` | `false` | not ported — see Fidelity (upstream default `true`) |

Parameter keys (`gatk_het_prior`, `deepvariant_model_type`,
`callable_sites_min_coverage`, `callable_sites_max_coverage`,
`callable_sites_fraction`, `callable_sites_merge_distance`,
`callable_sites_kmer`, `callable_sites_min_score`, `postprocess_contig_size`,
`postprocess_maf`, `postprocess_missingness`, `postprocess_exclude_scaffolds`,
`qc_min_depth`, `qc_max_sample_missingness`, `qc_exclude_scaffolds`,
`qc_clusters`, `qc_google_api_key`) mirror upstream's config values; pass
overrides as `oxo-flow run main.oxoflow key=value`.

## Source

Ported from **[harvardinformatics/snparcher](https://github.com/harvardinformatics/snparcher)**,
version `v2.2` (commit `e0e7a9478d4e042fce217db4e6077dafdaf57245`, MIT).
Created 2026-08-15; this workflow may lag behind upstream releases. Upstream
license and attribution are recorded in [NOTICE.md](NOTICE.md).

## Fidelity

58 rules ported from upstream v2.2 (up from 12), covering every branch that
can be expressed in oxo-flow's static DAG. Commands are ported verbatim
(same flags, same output paths); upstream's snakemake `{{...}}` shell escaping
is unwrapped to literal braces, and snakemake `{params.*}`/`{resources.*}`
references are resolved to their upstream values or to `{config.*}` keys.

| Upstream rule | oxo-flow rule | Tool (version) | Notes |
|---|---|---|---|
| `prepare_reference` (local branch) | `prepare_reference` | samtools 1.24 (bgzip) | identical command; url/accession branches not ported (config `reference_source` is a local path) |
| `index_reference` | `index_reference` | samtools 1.24, bwa 0.7.19 | identical command (faidx + dict + bwa index) |
| `fastp` | `fastp` / `fastp_srr` | fastp 1.3.6 | identical flags; `{sample}/{sample}/u1` fan-out for single-row sheets with empty `library_id`; SRA-downloaded reads via `fastp_srr` |
| `download_sra` | `download_sra` | sra-tools 3.2.1, ffq, curl, pigz | identical prefetch→fasterq-dump→pigz flow with ffq/ENA fallback; fasterq-dump `--tmpdir` dropped (oxo-flow has no per-rule tmpdir) |
| `bwa_mem` | `bwa_mem` | bwa 0.7.19, samtools 1.24 | identical command incl. read group `ID:{sample}.u1 SM:{sample} LB:{sample} PL:ILLUMINA`; raw BAM is temp like upstream |
| `merge_library_bams` | `merge_library_bams` | samtools 1.24 | per-library merge, single input unit in the default path |
| `merge_library_level_bams` | `merge_library_level_bams` | samtools 1.24 | no-markdup path (`results/bams/merged/{sample}.bam`) |
| `markdup_library` / `merge_dedup_libraries` | `markdup_library` / `merge_dedup_libraries` | sambamba 1.0.1, samtools 1.24 | identical commands; gated on `mark_duplicates` (default `false`; upstream default `true` — see deviations below) |
| `index_bam_csi` | `index_bam_csi` / `index_bam_csi_markdup` / `index_bam_csi_external` | samtools 1.24 | identical (`samtools index -c`), one per BAM-producing branch |
| `stage_external_bam` | `stage_external_bam` | — | external BAM inputs symlinked into `results/bams/input/` then run through the standard callers |
| `gatk_haplotypecaller` (standard mode) | `gatk_haplotypecaller` / `_markdup` / `_external` | gatk4 4.6.2.0 | identical flags incl. `-ploidy 2 --emit-ref-confidence GVCF --min-pruning 1 --min-dangling-branch-length 1` (low-coverage defaults); `-Xmx7000m` = upstream default profile `mem_mb_reduced`; threads 1 as upstream |
| `deepvariant_call` | `deepvariant_call` / `_markdup` / `_external` | google/deepvariant:1.10.0 (docker) | identical `/opt/deepvariant/bin/run_deepvariant` invocation; gated on `variant_tool = "deepvariant"` |
| `create_db_mapfile` | `create_db_mapfile` | python (script) | identical logic, ported as `scripts/write_joint_gvcf_mapfile.py` |
| `joint_genomics_db_import` | `joint_genomics_db_import` | gatk4 | identical GenomicsDBImport flow incl. `TILEDB_DISABLE_FILE_LOCKING` and `--merge-input-intervals` from `scripts/interval_list_tools.py` (merge threshold 50 = upstream `GENOMICSDB_MERGE_CONTIG_THRESHOLD`) |
| `joint_genotype_gvcfs` | `joint_genotype_gvcfs` | gatk4 | identical (tar-extract → `gendb://` GenotypeGVCFs → `results/vcfs/raw.vcf.gz`); temp raw VCF like upstream |
| `glnexus_joint` | `glnexus_joint` | glnexus, bcftools 1.23 | identical DeepVariant-config GLnexus join; `mem_gbytes` = `mem_mb_reduced/1024` rounded to 8, computed from the default profile (see deviations) |
| `variant_filtration` | `variant_filtration` | gatk4, bcftools | identical RPRS/FS_SOR/MQ/QUAL hard filters, `--invalidate-previous-filters true`, then `bcftools index -f -t` |
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
| qc module (`contig_map`, `vcftools_individuals`, `subsample_snps`, `prepare_plink_inputs`, `copy_qc_report`, `plink`, `setup_admixture`, `admixture`, `qc_dashboard`) | `qc_contig_map` … `qc_dashboard` | vcftools 0.1.16, bcftools 1.23, plink2/plink, admixture, R | identical commands; logic ported to `scripts/contig_map.py`, `vcftools_individuals.py`, `prepare_plink_inputs.py`, `contigs4admixture.py`, `qc_dashboard_render.R` |
| `setup` / `download_reads` / `map_samples` / `call_variants` / `qc_report` / `callable_sites` / `gvcfs` (Snakefile aggregation targets) | n/a | — | Snakemake target rules, no commands of their own |

### Remaining exclusions (structurally impossible in oxo-flow)

| Item | Why excluded | Evidence |
|---|---|---|
| `intervals.enabled` interval scatter (`picard_intervals`, `create_gvcf_intervals`, `create_db_intervals`, `gatk_haplotypecaller_interval`, `concat_interval_gvcfs*`, `concat_interval_vcfs*`, `compress_interval_raw_vcf`, `normalize_external_gvcf_for_gatk`, `archive_gatk_gvcf`) | upstream uses snakemake **checkpoints** (`create_gvcf_intervals` at `workflow/rules/intervals.smk:59`, `create_db_intervals` at `:89`) that extend the DAG at runtime from intermediate files; oxo-flow's DAG is static, so per-interval fan-out cannot be planned | `intervals.smk:59,89` (`checkpoint`), `variant_calling/gatk.smk` consumes `intervals.enabled` via checkpoints |
| `bcftools_call` (bcftools caller) | depends on `bcftools_regions` **checkpoint** (`variant_calling/bcftools.smk:13`) which computes regions from the runtime `callable_sites.bed`; same static-DAG limitation, and the checkpoint also writes the per-sample gVCF into the same output paths the GATK/DeepVariant branches produce (runtime producer selection) | `variant_calling/bcftools.smk:13` (`checkpoint bcftools_regions`), `:42,:55` |
| parabricks (all `parabricks_*` rules) | requires `--nv` GPU passthrough (upstream `parabricks.smk` runs `--nv` images with `nvidia-docker`); the oxo-flow docker backend has no `--nv` support and no GPU device declaration; additionally NVIDIA EULA/license enforcement cannot be guaranteed in CI | `variant_calling/parabricks.smk` (every rule is `--nv`) |
| sentieon (all `sentieon_*` rules) | proprietary tool gated on a `SENTIEON_LICENSE` server and a pre-installed license; cannot be distributed or verified in a community port | `config/config.yaml` `sentieon` section; `workflow/rules/sentieon.smk` |
| `denovo` and `structural_variants` pipeline sections | do not exist as rules in upstream v2.2 | grep of `workflow/` at e0e7a94 finds neither rule set |
| `gvcf` input type (`normalize_external_gvcf_for_gatk`) | external per-sample paths (arbitrary filesystem locations) cannot feed `expand_inputs` group expansion, which requires a uniform `{sample}` pattern under the workflow root | `variant_calling/gatk.smk:57` (`rule normalize_external_gvcf_for_gatk`) |
| `generate_coords_file` (qc module) | upstream generates a `coords.tsv` from sample sheet `lat`/`long` columns; the port's sample-group model has no lat/long metadata, and it only feeds the excluded gvcf-input dashboard variant | `modules/qc/Snakefile` `generate_coords_file` |
| per-sample `mark_duplicates` override | upstream reads the value from the sample sheet per row; the port maps it to the global `mark_duplicates` config key | `config/config.yaml` `mark_duplicates` vs `workflow/snakefiles` per-sample handling |
| multi-library / multi-unit rows (library_id, input_unit) | the sample-group model has one unit per sample; consumers of `results/bams/raw/{sample}/{sample}/u1.bam` are hard-coded to the `u1` unit | sample sheet semantics in upstream `README` |

### Documented deviations from upstream

1. **Default config differs from upstream**: upstream defaults `mark_duplicates: true`, `joint_genotyping: true`, `generate_filtered_vcf: true`, `callable_sites.enabled: true`; the port defaults all of these to `false` so its default plan is byte-identical to the previous 12-rule port. Flip the keys to get upstream's full pipeline.
2. **postprocess/qc modules consume `results/vcfs/raw.vcf.gz`**, not upstream's `FINAL_VCF` (which is the hard-filtered VCF when `generate_filtered_vcf: true` and GATK is used). The modules were run on the raw joint VCF. The difference only matters when combining GATK + `generate_filtered_vcf` + a module, and reproduces upstream behavior for the DeepVariant path.
3. **Long-contig CSI mode not ported for postprocess**: upstream conditionally uses CSI indexes when contigs exceed 512 Mb (`regions_to_index`); the port always uses the default TBI short mode. Applicable only to genomes with >512 Mb contigs.
4. **`glnexus_joint` memory**: upstream computes `mem_gbytes` from the default profile's `mem_mb_reduced`; the port inlines the resulting value 8 (with the same `if < 1 then 1` clamp).
5. **`combine_qc_metrics` with mixed input types**: its `expand_inputs` references `results/fastp/{sample}/{sample}/u1.json` for every sample, which only exists for fastq/srr samples — for bam-input cohorts the fastp stats expand fails at plan time. Use fastq/srr groups when combining QC metrics.
6. **`fasterq-dump --tmpdir` dropped** (no per-rule tmpdir in oxo-flow); SRA downloads use the current directory.

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
