# Percolator vendored source tree

- **Upstream**: https://github.com/percolator/percolator
- **Pinned at**: commit eb157f7 (post-rel-3-08-01; PR #399 / I-spline PEP)
- **Commit SHA**: eb157f74e963430e559e0d0bcd31291e4ad660ba
- **License**: Apache 2.0 (see LICENSE-Apache-2.0.txt) plus the embedded SVMlin files (see NOTICE-percolator.txt).

## What's here

A stripped-down copy of Percolator's `src/` tree, covering the PSM rescoring path only:
XML / tab / CLI / protein inference / peptide fragmentation code is excluded.

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

## Optional: 3.08 parity tests

`Percolator_subprocess_parity_test`'s §9–§10 sections compare in-process
`pep_method="isotonic"` / `"logistic_regression"` against a 3.08-era
subprocess `percolator` binary (which exposes `--pava-pep` / `--ip-pep` /
`--irls-pep`; the bundled 3.06 binary predates these flags). The sections
skip silently unless a 3.08 binary is available.

To enable them locally:

```bash
# Fetch the rel-3-08 prebuilt (no-XML) deb and extract just the binary
mkdir -p OpenMS-build/THIRDPARTY/Percolator-eb157f7
curl -sL -o /tmp/p38.deb \
  https://github.com/percolator/percolator/releases/download/rel-3-08/percolator-noxml-v3-08-linux-amd64.deb
dpkg-deb -x /tmp/p38.deb /tmp/p38-extracted
cp /tmp/p38-extracted/usr/bin/percolator OpenMS-build/THIRDPARTY/Percolator-eb157f7/
cmake --build OpenMS-build -j$(nproc) --target OpenMS  # re-runs configure, picks up the binary
```

CMake prints `Percolator_subprocess_parity_test: 3.08 parity sections
enabled` on configure when the binary is found.

Note: rel-3-08 (3.08.0, May 2025) is slightly older than the eb157f7 SHA
we vendor, so observed drift is ~0.28 on PEP/q-value even under
same-algorithm comparison. This is binary-version skew, not a parity
regression; test tolerances (0.5) accommodate it. The `r ≥ 0.9999`
assertion on SVM scores is the strong signal — any regression that
corrupts the SVM path fails it immediately.

## Licensing note

The file-level license of the vendored `ssl.{h,cpp}` (SVMlin) files is
understood to have been explicitly resolved to Apache-2.0-compatible terms by
the OpenMS maintainers prior to landing this subtree. See `NOTICE-percolator.txt`
for the attribution.
