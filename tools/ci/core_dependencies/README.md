# Proposed core dependency boundary

This implements item 1 of [#10083](https://github.com/OpenMS/OpenMS/issues/10083).
It provides an include inventory and prevents new crossings of a proposed core
boundary while the existing couplings are removed in subsequent PRs.

Run from the repository root with Python 3.9 or newer. No third-party Python
packages, CMake configuration, compiler, or OpenMS build are required:

```sh
python -m unittest discover -s tools/ci -p 'test_check_core_dependencies.py' -v
python tools/ci/check_core_dependencies.py --report core-dependencies.json
python tools/ci/check_core_dependencies.py --update-core   # after changing the boundary
```

The GitHub Actions workflow runs on every PR targeting `develop`, pushes to
`develop`, and manual dispatch. The JSON report is uploaded even if the boundary
check fails. There are no path filters, so this lightweight job can be used as a
required check without unrelated PRs leaving it pending.

## Boundary and exceptions

The core is derived, not classified by hand. Starting from five installed
headers (`Peak1D`, `ChromatogramPeak`, `MSSpectrum`, `MSChromatogram`,
`MSExperiment`), the scanner follows header includes, each header's mirrored
implementation file (`include/OpenMS/X/A.h` -> `source/X/A.cpp`), and that
file's includes. It does not follow the edges recorded in `baseline.json`; those
are the boundary. The result is written to `core_files.txt` by `--update-core`
and checked in, so the check can name each crossing edge instead of silently
absorbing whatever a new include pulls in.

The first snapshot has 67 headers and 59 implementation files. This is a
candidate boundary, not a declaration that all these files already form an
independently linkable library. In particular, `File.cpp`, `MSSpectrum.cpp`, and
`DocumentIdentifier.cpp` still contain operations requiring higher-level services.

The check enforces these rules against the snapshot:

| Rule | Forbidden edge |
| --- | --- |
| `core-to-noncore` | A core file includes a resolved project file outside the core. |
| `core-unresolved` | A core file has an unresolved project/quoted include or a macro include that cannot be classified. |
| `installed-to-private` | Any installed header in the scanned libraries includes a resolved, non-installed project file. |

`baseline.json` records existing violations by **rule, source file, and target
file**, with a reason for each exception. It starts with ten edges. Line numbers
and duplicate include occurrences do not change an exception's identity. There
are no directory exemptions and no command that silently accepts new violations.
New violations and stale baseline entries both fail the check, as do differences
between `core_files.txt` and the derived closure. Configuration or checkout
errors return exit code 2; boundary/baseline/snapshot failures return 1.

Workflow for a failing check:

- A new `core-to-noncore` edge: remove the include, or add a baseline entry with
  a reason. Either way the derived closure no longer changes.
- The include is a legitimate new core helper: run `--update-core`, review the
  files it prints, and commit the regenerated snapshot. If the helper's own
  implementation reaches outside the boundary, baseline those edges first so the
  snapshot grows by the helper alone.
- A removed dependency: delete its exact baseline entry in the same PR and run
  `--update-core` if the closure shrank.

Only the implementation file mirroring a core header is part of the core. When
methods are split across translation units, the split-off file (for example a
future `FileConfig.cpp` next to `File.cpp`) stays outside the boundary, and its
includes are not checked. Reviewers must check that the methods left in the
mirrored file are the ones meant to be core: textual scanning cannot identify
which methods a translation unit defines.

## Inventory

The scanner reads headers and implementation files under `include/` and `source/`
in `src/openms` and `src/openswathalgo`. Vendored dependencies, GUI, TOPP, bindings,
and tests are outside this initial inventory. The deterministic JSON contains:

- Per-file kind, installation status, and core membership.
- Separate `header_includes` and `implementation_includes`, with original spelling,
  source line, resolved target, and status. `OpenMS/...`, `include/OpenMS/...`, and
  relative quoted includes are recognized. All conditional branches are included,
  including `#if 0`. Comments and raw-string examples are excluded.
- `internal`, `generated`, `external`, `unresolved`, and `macro` include statuses.
  Only the generated headers listed in the script receive that status.
- Reachability from each seed header, following include edges only, and the
  derived `core_files` closure, which additionally steps from a header to its
  mirrored implementation file.
- File-level strongly connected components in `cycles`. An empty list means no
  textual include cycles were found; it says nothing about symbol/link cycles.
- A separate `runtime_resources` inventory of literal `File::find("...")`
  arguments. This includes the PSI-MS, PATO, UO, BTO, and GO vocabulary files loaded
  by `ControlledVocabulary.cpp`, even though that loader is outside the core.
- Current violations, new violations, stale baseline entries, and snapshot
  differences (`new_core_files`, `stale_core_files`).

Installed headers come from the `sources_list_h` lists registered in
`src/openms/includes.cmake` and the `header_algo_list`/`header_dataaccess_list`
lists in `OpenSwathAlgoFiles.cmake`, which feed the libraries' `HEADER_FILES`
arguments. Conditional lists are unioned across configurations; explicitly
private or unlisted headers are not marked installed. Generated headers are
reported separately because their output files need not exist in a source checkout.

## Limits

This is source scanning, not preprocessing or symbol analysis. It does not prove
buildability, link independence, ABI compatibility, or unused includes. It does
not expand include macros, evaluate feature flags, search compiler/system include
paths, or parse arbitrary CMake programs. Unresolved quoted includes may refer to
external dependencies supplied by the build. System and third-party angle includes
are inventoried as external; external-library policy is not enforced in this first
guard. A future external-dependency rule should distinguish standard/platform
headers from libraries rather than guessing from include spelling.

The CMake reader supports literal `set` and `list(APPEND)` values in the current
install-list conventions. It rejects computed list values and unsupported list
mutations. Update the reader and its tests if those conventions change.

Runtime data discovery is intentionally partial: it records literal first
arguments, including path prefixes, but cannot resolve dynamically constructed
filenames, resources loaded through other APIs, or the runtime call graph. These
entries must not be interpreted as a complete deployment manifest.
