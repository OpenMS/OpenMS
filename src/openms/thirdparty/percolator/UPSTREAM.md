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
  drops, etc. Regenerated 2026-04-24 (commit 9f227d7c03) from the
  current-tree-vs-upstream diff — replaces the former 01-06 chain which
  had accumulated drift.

## Licensing note

The file-level license of the vendored `ssl.{h,cpp}` (SVMlin) files is
understood to have been explicitly resolved to Apache-2.0-compatible terms by
the OpenMS maintainers prior to landing this subtree. See `NOTICE-percolator.txt`
for the attribution.
