# Complete libxml2 replacement for OpenMS 3.7

This feature branch starts at develop `d7d9974a63` (OpenMS 3.7.0) and carries the complete native libxml2 migration, including the reviewed optimizations from `ac0e80acbe`. It deliberately leaves the prototype branches and benchmark archives in their original history. The upstream chromatogram-skipping fix and current HiGHS solver/dependency defaults are retained.

## Implementation

- Replace Xerces SAX infrastructure, attributes/string conversion, XML buffers, compressed input streams, schema validation, mzIdentML DOM handling and imzML parsing. Public headers no longer include Xerces types; downstream subclasses of the old parser-specific interfaces may need source changes.
- Require libxml2 2.13 or newer with SAX1, PUSH, READER, SCHEMAS and THREAD support. Discovery checks those features and warns if both iconv and ICU are absent. Dependencies use `LibXml2::LibXml2`. Discovery prefers libxml2's CMake configuration, which carries static compile definitions and private dependencies, then falls back to CMake's Find module for installations without that configuration. OpenMS itself remains a shared library. Custom static libxml2 installations should provide their upstream configuration or equivalent dependency metadata.
- Preserve the reviewed namespace, character/encoding, entity, early-stop and diagnostic behavior. External entities remain disabled. The documented schema-instance compatibility difference remains: libxml2 accepts an invalid-URI fixture that Xerces rejects, although all bundled schemas compile. This is not a claim of identical validator behavior for every document.
- Keep native UTF-8 handling for opted-in handlers, avoid transcoding skipped binary payloads, compare ASCII attribute names without conversion, and use checked binary reads for indexed footers.
- For wrapped index footers above 2 MiB, traverse forward and expand only an individual offset subtree. Smaller footers retain the DOM path to avoid the measured ion-mobility read-loop regression. Both paths preserve index semantics, transactional failure outputs and aliased-output ordering. No process-wide allocator tuning or trim call is introduced. On ELF, static libxml2 symbols stay private to OpenMS to avoid interposition with other XML runtimes loaded by the application.

## Build and distribution support

The vcpkg manifest removes `xerces-c` and raises the existing libxml2 dependency minimum. The current vcpkg baseline/submodule and overlay dependency sources are unchanged. The shared CI build action disables Xerces discovery and checks source/manifests/native-library dependencies after building. The existing Build and Test workflow covers Windows MSVC x64, macOS Apple Clang ARM64, and Linux GCC/Clang x64 and GCC ARM64. The wheel workflow and Docker dependency recipe use libxml2; the Docker source download is checksum-pinned because the base distribution's version is older than the required minimum.

Use the existing platform presets, for example:

```sh
git submodule update --init vcpkg
cmake --preset linux-x64-release
cmake --build --preset linux-x64-release
```

Select the corresponding macOS or Windows preset for those systems. System-package builds can set `OPENMS_USE_VCPKG=OFF`, with a sufficiently recent libxml2 prefix. Windows builds require the supported Visual Studio toolchain; static dependency builds must match the selected MSVC runtime/configuration. GUI, wheels and containers retain their existing platform requirements. With vcpkg, the deprecated COIN/GLPK solver options additionally need the `lp-coin`/`lp-glpk` manifest feature.

Regression guards:

```sh
python -m unittest discover -s tools/ci/tests -p test_check_no_xerces.py
python tools/ci/check_no_xerces.py --build-dir build/linux-x64-release
```

The guard inspects built OpenMS libraries using the native platform tool, and fails if inspection is unavailable or a Xerces runtime dependency is present. It does not promise that unrelated optional third-party programs never use another XML library.

## Performance and prior validation

The [retained-context report and graphs](https://github.com/okohlbacher/OpenMS/blob/4d28743d348ef893a00c6049f23d3d4d2b75c0d6/doc/xml-libxml2-validation/context-memory/REPORT.md) preserve all measured results, identities, reviews and rejected iterations. At 96 Linux readers the selected strategy reduced retained index/context memory by 72–87% on four inputs, left the small ion-mobility input unchanged, and had no sustained read-loop median regression above the 5% investigation gate (worst +1.8%). Teardown and post-close RSS can still be worse; the report records those costs separately. Those benchmarks belong to their original product source and were not repeated merely to integrate newer develop changes.

The [earlier complete-migration evidence](https://github.com/okohlbacher/OpenMS/blob/4d28743d348ef893a00c6049f23d3d4d2b75c0d6/doc/xml-libxml2-validation/README.md) contains Linux/macOS full-suite results, schema/API/encoding controls and Codex/Kimi/Claude Opus reviews. They are historical evidence, not a successful Windows or current-base build claim. Fresh build/test outcomes for this branch are recorded below when completed; workflow configuration alone is not a passing build.
