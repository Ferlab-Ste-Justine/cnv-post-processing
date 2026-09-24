#!/usr/bin/env bash
# Automated pre-push checks: pre-commit, nf-test suite, and nf-core lint. Fail-fast
# (set -e) so the first broken step stops the script.
#
# pre-commit (prettier, trailing-whitespace, end-of-file-fixer, config in
# .pre-commit-config.yaml) also runs in CI on every PR
# (.github/workflows/linting.yml), but it's included here too so it fails
# fast locally instead of only on push.
#
# Note on dependencies: the dry-run and lint steps need nothing but the
# tools themselves, but the real (non-dry-run) nf-test suite still spins up
# Docker containers, and its pipeline-level tests (tests/default.nf.test)
# need the S3 test dataset synced locally (data-test/ -- see CLAUDE.md's
# "Test dataset" section). This script does NOT include the ad hoc, manually
# -inspected `nextflow run` smoke tests -- see scripts/run-smoke-tests.sh for
# those.
#
# Usage: scripts/run-test-suite.sh

set -euo pipefail

cd "$(git rev-parse --show-toplevel)"
export NXF_FILE_ROOT="$PWD"

require() {
    command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found on PATH." >&2; exit 1; }
}
require pre-commit
require nf-test
require nf-core
require node

step() { echo; echo "==> $1"; }

# If .snap file needs to be generated.
# nf-test test modules/local/slivar_expr --updateSnapshot

step "[1/5] pre-commit (prettier, trailing-whitespace, end-of-file-fixer)"
pre-commit run --all-files

step "[2/5] slivar classification function unit tests (plain node, no Docker needed)"
node assets/tests/cnv_slivar_functions.test.js

step "[3/5] nf-test dry-run (syntax check only, no execution)"
nf-test test --dryRun

step "[4/5] nf-test full suite (--profile test,docker)"
nf-test test --profile test,docker

step "[5/5] nf-core pipelines lint --release"
nf-core pipelines lint --release

echo
echo "All checks passed."
