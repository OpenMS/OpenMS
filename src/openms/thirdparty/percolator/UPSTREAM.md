# Percolator vendored source tree

- **Upstream**: https://github.com/percolator/percolator
- **Pinned at**: rel-3-06-1 (update via sync-from-upstream.sh)
- **Commit SHA**: (will be filled by sync-from-upstream.sh)
- **License**: Apache 2.0 (see LICENSE-Apache-2.0.txt) plus the embedded SVMlin files (see NOTICE-percolator.txt).

## What's here

A stripped-down copy of Percolator's `src/` tree, covering the PSM rescoring path only:
XML / tab / CLI / protein inference / peptide fragmentation code is excluded.

## How to re-sync

Run `./sync-from-upstream.sh <tag-or-sha>`. The script copies files listed in
`whitelist.txt`, applies `patches/*.patch`, and records the new upstream commit
SHA in this file. Review the resulting diff and commit as
`chore(percolator): sync from upstream <ref>`.

## Licensing note

The file-level license of the vendored `ssl.{h,cpp}` (SVMlin) files is
understood to have been explicitly resolved to Apache-2.0-compatible terms by
the OpenMS maintainers prior to landing this subtree. See `NOTICE-percolator.txt`
for the attribution.
