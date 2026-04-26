# Percolator vendored source tree

- **Upstream**: https://github.com/percolator/percolator
- **Pinned at**: commit eb157f7 (post-rel-3-08-01; PR #399 / I-spline PEP)
- **Commit SHA**: eb157f74e963430e559e0d0bcd31291e4ad660ba
- **License**: Apache 2.0 (see LICENSE-Apache-2.0.txt) plus the BSD 3-Clause licensed liblinear-derived TRON solver (see NOTICE-percolator.txt).

## What's here

A stripped-down copy of Percolator's `src/` tree, covering the PSM rescoring path only:
XML / tab / CLI / protein inference / peptide fragmentation code is excluded.

### Local replacement: TRON-based SVM solver

`ssl.{cpp,h}` and `tron.{cpp,h}` come from the
[`percolator-tron`](https://bitbucket.org/jthalloran/percolator_upgrade)
development branch (Halloran's TRON integration of liblinear v2.11). They
replace the legacy SVMlin-based L2-SVM-MFN solver from upstream
Percolator's `ssl.{cpp,h}`. The replacement:
- Sidesteps the SVMlin license question (its upstream headline is GPL v2+;
  the relicense to Apache-2.0 was scoped to Percolator itself, not to
  link-time inclusion in downstream non-Apache projects).
- Calls `extern "C"` BLAS (`dnrm2_`, `ddot_`, `daxpy_`, `dscal_`) which
  resolve at link time against OpenMS's existing libblas dependency.
- Is NOT in `whitelist.txt` — `sync-from-upstream.sh` won't overwrite
  these files. If a future upstream Percolator sync introduces a different
  ssl.cpp, this replacement will need to be re-applied manually.

## How to re-sync

Run `./sync-from-upstream.sh <tag-or-sha>`. The script copies files listed in
`whitelist.txt`, applies `patches/*.patch`, and records the new upstream commit
SHA in this file. Review the resulting diff and commit as
`chore(percolator): sync from upstream <ref>`.

A clean re-sync of the currently-pinned SHA reproduces the committed tree
bit-identically (verified 2026-04-24). If a re-sync produces a non-empty diff
(beyond `UPSTREAM_COMMIT` timestamp churn), it means either the pinned SHA
moved or a new downstream adaptation is needed — regenerate
`patches/01-namespace-wrap.patch` from the post-adaptation tree rather than
hand-editing sources and committing them, to keep sync self-healing.

## Patches

- `patches/01-namespace-wrap.patch` — single comprehensive patch covering all
  OpenMS adaptations: namespace wrap in `OpenMS::Internal::Percolator`,
  `using namespace std;` placement, `std::` qualifications, additional
  includes (TabReader, Version.h), `::min/::max → std::min/std::max`,
  `parseOptions` body removal, enzyme / protein-inference / SQT feature
  drops, `#pragma once` on the new header-only regressors, etc.
  Regenerated 2026-04-24 (commit 9f227d7c03) from the current-tree-vs-upstream
  diff — replaces the former 01-06 chain which had accumulated drift.

### Regenerating the patch

Whenever you hand-fix a downstream issue in the vendored tree, capture it in
this patch rather than leaving it only in the source commit. To regenerate:

```bash
mkdir -p /tmp/perc-eb157f7
whitelist=$(grep -v '^#' whitelist.txt | grep -v '^$')
for f in $whitelist; do
  gh api "repos/percolator/percolator/contents/src/$f?ref=$(cat UPSTREAM_COMMIT)" \
    --jq '.content' 2>/dev/null | base64 -d > "/tmp/perc-eb157f7/$f" 2>/dev/null
done
: > patches/01-namespace-wrap.patch
for f in $(ls /tmp/perc-eb157f7/ | sort); do
  [ -s "/tmp/perc-eb157f7/$f" ] || continue
  diff -urN --label "/tmp/percolator-preserved/$f" --label "src/openms/thirdparty/percolator/$f" \
    "/tmp/perc-eb157f7/$f" "src/openms/thirdparty/percolator/$f" \
    >> patches/01-namespace-wrap.patch || true
done
```

The `--label` flags are important — without them, `diff -urN` embeds
filesystem timestamps that make `git diff` noisy on every regeneration even
when the semantic patch content is unchanged.

## Licensing note

`ssl.{cpp,h}` and `tron.{cpp,h}` are BSD-3-Clause (liblinear v2.11 via
the percolator-tron fork). All other vendored files are Apache-2.0
(Percolator). See `NOTICE-percolator.txt` for the full attribution.
