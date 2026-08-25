oxo-flow-snparcher
Copyright (c) 2026 oxo-flow-community

This pipeline is a port of harvardinformatics/snparcher
(https://github.com/harvardinformatics/snparcher), version v2.2
(commit e0e7a9478d4e042fce217db4e6077dafdaf57245), authored by
Harvard Informatics.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

---------------------------------------------------------------------
Upstream license

This port is derived from harvardinformatics/snparcher under the MIT
license. The upstream LICENSE is included verbatim in this repository
at LICENSE.upstream (fetched from the upstream repository at the ported
commit).
---------------------------------------------------------------------
Script attribution

The following files are copied verbatim from the upstream repository at
the ported commit (diff-verified byte-identical); they retain their
upstream authorship and are used under the upstream MIT license:

- scripts/callable_coverage_thresholds.py
  (upstream: workflow/scripts/callable_coverage_thresholds.py)
- scripts/callable_zarr_to_bed.py
  (upstream: workflow/scripts/callable_zarr_to_bed.py)
- scripts/qc_dashboard_interactive.Rmd
  (upstream: workflow/modules/qc/scripts/qc_dashboard_interactive.Rmd)
- scripts/interval_list_tools.py
  (upstream: workflow/scripts/interval_list_tools.py)

The following files are adapted from upstream scripts (same logic,
repackaged for the oxo-flow execution model — CLI arguments, file
handling, and subprocess plumbing differ; logic is otherwise identical):

- scripts/contig_map.py — plain-stdlib port of the upstream
  workflow/modules/qc/Snakefile `contig_map` rule's pandas run block
- scripts/contigs4admixture.py — CLI adaptation of the upstream
  workflow/modules/qc/Snakefile `setup_admixture` run block
- scripts/qc_dashboard_render.R — CLI adaptation of the upstream
  workflow/modules/qc/Snakefile `qc_dashboard` run block (calls
  rmarkdown::render on the verbatim qc_dashboard_interactive.Rmd)
- scripts/vcftools_individuals.py — CLI adaptation of the upstream
  workflow/modules/qc/Snakefile `vcftools_individuals` run block
- scripts/prepare_plink_inputs.py — CLI adaptation of the upstream
  workflow/modules/qc/Snakefile `prepare_plink_inputs` run block

All scripts retain upstream authorship and are used under the upstream
MIT license.
---------------------------------------------------------------------
