# Proposed core dependency boundary

This implements item 1 of [#10083](https://github.com/OpenMS/OpenMS/issues/10083).
It provides an include inventory and prevents new crossings of a proposed core
boundary while the existing couplings are removed in subsequent PRs.

Run from the repository root with Python 3.9 or newer. No third-party Python
packages, CMake configuration, compiler, or OpenMS build are required:

```sh
python -m unittest discover -s tools/ci -p 'test_check_core_dependencies.py' -v
python tools/ci/check_core_dependencies.py --report core-dependencies.json
```

The GitHub Actions workflow runs on every PR targeting `develop`, pushes to
`develop`, and manual dispatch. The JSON report is uploaded even if the boundary
check fails. There are no path filters, so this lightweight job can be used as a
required check without unrelated PRs leaving it pending.

## Boundary and exceptions

`manifest.json` starts from five installed headers: `Peak1D`, `ChromatogramPeak`,
`MSSpectrum`, `MSChromatogram`, and `MSExperiment`. It explicitly classifies each
reachable header helper and the implementation files selected for the candidate
core. Reasons explain inclusion of value types, metadata, logging, path utilities,
file-type identifiers, and ion-mobility identifiers. For example, `FileTypes` is
inside the boundary and `FileHandler` is outside it. Directory names alone do not
decide whether a file belongs to the core.

The first manifest has 67 headers and 59 implementation files. This is a candidate
boundary, not a declaration that all these files already form an independently
linkable library. In particular, `File.cpp`, `MSSpectrum.cpp`, and
`DocumentIdentifier.cpp` still contain operations requiring higher-level services.

The check enforces these rules:

| Rule | Forbidden edge |
| --- | --- |
| `core-to-noncore` | A classified core file includes a resolved project file outside the core. |
| `core-unresolved` | A classified core file has an unresolved project/quoted include or a macro include that cannot be classified. |
| `installed-to-private` | Any installed header in the scanned libraries includes a resolved, non-installed project file. |

`baseline.json` records existing violations by **rule, source file, and target
file**, with a reason for each exception. It starts with ten edges. Line numbers
and duplicate include occurrences do not change an exception's identity. There
are no directory exemptions and no command that silently accepts new violations.
New violations and stale baseline entries both fail the check. Configuration or
checkout errors return exit code 2; boundary/baseline failures return 1.

When removing a dependency, remove its exact baseline entry in the same PR. When
introducing a legitimate core helper, review its dependencies and classify the
header and each relevant implementation file explicitly, with a reason. A new
non-core helper does not inherit an exception from its directory or class name.

When splitting methods across translation units, classify each resulting file
according to its role. The scanner never combines `File.h`, `File.cpp`, and a
future `FileConfig.cpp` into one class node. Keep basic path operations inside the
candidate core; a configuration implementation can remain outside it. Reviewers
must check this classification: textual scanning cannot identify which methods a
new implementation file defines.

## Inventory

The scanner reads headers and implementation files under `include/` and `source/`
in `src/openms` and `src/openswathalgo`. Vendored dependencies, GUI, TOPP, bindings,
and tests are outside this initial inventory. The deterministic JSON contains:

- Per-file kind, installation status, and core classification reason.
- Separate `header_includes` and `implementation_includes`, with original spelling,
  source line, resolved target, and status. `OpenMS/...`, `include/OpenMS/...`, and
  relative quoted includes are recognized. All conditional branches are included,
  including `#if 0`. Comments and raw-string examples are excluded.
- `internal`, `generated`, `external`, `unresolved`, and `macro` include statuses.
  Only the exact generated headers listed in the manifest receive that status.
- Reachability from each seed header and from all explicitly classified core
  files. These walks follow include edges only; they never invent a header-to-cpp
  edge. A header's include closure therefore stays separate from its implementation.
- File-level strongly connected components in `cycles`. An empty list means no
  textual include cycles were found; it says nothing about symbol/link cycles.
- A separate `runtime_resources` inventory of literal `File::find("...")`
  arguments. This includes the PSI-MS, PATO, UO, BTO, and GO vocabulary files loaded
  by `ControlledVocabulary.cpp`, even though that loader is outside the core.
- Current violations, new violations, and stale baseline entries.

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
