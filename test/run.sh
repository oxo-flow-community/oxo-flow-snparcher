#!/usr/bin/env bash
# Acceptance test for the oxo-flow-snparcher port.
# Usage: ./test/run.sh            (uses ./main.oxoflow)
set -euo pipefail
cd "$(dirname "$0")/.."
OXO=${OXO:-oxo-flow}

echo "==> validate"
"$OXO" validate main.oxoflow

echo "==> lint (warnings are acceptable, errors are not)"
"$OXO" lint main.oxoflow

echo "==> dry-run with default config"
# oxo-flow v0.11.0 prints the plan to stderr; capture both streams
"$OXO" dry-run main.oxoflow --samples first:1 > /tmp/oxo-dryrun-$$.txt 2>&1
grep -q "would execute" /tmp/oxo-dryrun-$$.txt

echo "==> per-sample mark_duplicates override dry-run (metadata_file)"
md=$(mktemp -d)
awk '/^author =/{print; print "metadata_file = \"samples.tsv\""; next} {print}' main.oxoflow > "$md/main.oxoflow"
printf 'sample\tmark_duplicates\nsample1\ttrue\nsample2\tfalse\n' > "$md/samples.tsv"
"$OXO" dry-run "$md/main.oxoflow" --samples first:2 mark_duplicates=false > /tmp/oxo-meta-$$.txt 2>&1
grep -q "markdup_library_cohort_sample1.*\[run" /tmp/oxo-meta-$$.txt || { echo "per-sample markdup override not applied (sample1)"; exit 1; }
grep -q "merge_library_level_bams_cohort_sample2.*\[run" /tmp/oxo-meta-$$.txt || { echo "per-sample nomarkdup override not applied (sample2)"; exit 1; }
# The row wins BOTH directions: with config mark_duplicates=true, sample2's
# row `false` must still take the no-markdup path (global true overridden).
"$OXO" dry-run "$md/main.oxoflow" --samples first:2 mark_duplicates=true > /tmp/oxo-meta-true-$$.txt 2>&1
grep -q "markdup_library_cohort_sample1.*\[run" /tmp/oxo-meta-true-$$.txt || { echo "per-sample markdup override not applied (sample1, config true)"; exit 1; }
grep -q "merge_library_level_bams_cohort_sample2.*\[run" /tmp/oxo-meta-true-$$.txt || { echo "row value must win over global true (sample2)"; exit 1; }
! grep -q "markdup_library_cohort_sample2.*\[run" /tmp/oxo-meta-true-$$.txt || { echo "row false must close the markdup gate despite global true (sample2)"; exit 1; }
# Any other value (e.g. "foo") must NOT override: it falls back to the global
# key — with config false, sample2 takes the no-markdup path and the markdup
# gate must NOT also open (no double-run).
printf 'sample\tmark_duplicates\nsample1\ttrue\nsample2\tfoo\n' > "$md/samples.tsv"
"$OXO" dry-run "$md/main.oxoflow" --samples first:2 > /tmp/oxo-meta-foo-$$.txt 2>&1
grep -q "markdup_library_cohort_sample1.*\[run" /tmp/oxo-meta-foo-$$.txt || { echo "per-sample markdup override not applied (sample1, non-boolean value case)"; exit 1; }
grep -q "merge_library_level_bams_cohort_sample2.*\[run" /tmp/oxo-meta-foo-$$.txt || { echo "non-boolean value must fall back to global (sample2)"; exit 1; }
! grep -q "markdup_library_cohort_sample2.*\[run" /tmp/oxo-meta-foo-$$.txt || { echo "non-boolean value must not open the markdup gate (sample2 double-run)"; exit 1; }
rm -rf "$md"

echo "==> debug: expanded commands contain no literal {wildcards}"
"$OXO" debug main.oxoflow 2>&1 | grep -q '{sample}' && { echo "unexpanded wildcards in debug output"; exit 1; } || true

echo "PASS"
