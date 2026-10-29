# Vendored mzPeak C++ library

Source: https://github.com/okohlbacher/mzpeak-openms, branch `trunk`, commit
`1d487a3a22be18d1e5fe9cf5f04517ceced7fd59` (2 Oct 2026), a fork of the upstream
https://github.com/OpenMS/mzpeak (`trunk`, `d70ce3ac21ff885efb23bb3b49f148907e7a5e3b`,
14 Sep 2026, the fork's merge base; the fork is 266 commits ahead, upstream has no
commits the fork lacks). License: BSD-3-Clause (`LICENSE`); the third-party notes of the
library are in `THIRD_PARTY.upstream.md`.

The library implements the mzPeak format of HUPO-PSI (https://github.com/HUPO-PSI/mzPeak):
`metadata.version` `0.9.0`, audited against `HUPO-PSI/mzPeak-specification` up to
commit `840ba4a` (29 Sep 2026) per the fork's `docs/spec-compliance.md`. The reader reads
both metadata layouts (nested and split facets), point and chunked signals (delta,
MS-Numpress), Bruker ims-compact, per-peak ion mobility; retention times are returned in
seconds. The writer emits the split-metadata layout.

Vendored: `include/` and `src/` of the library (not its tests, tools, examples, docs or
meson files). Built by `CMakeLists.txt` here as the static library `OpenMS_mzpeak`, linked
PRIVATE into libOpenMS when `OPENMS_ENABLE_MZPEAK` is ON; `OpenMS/FORMAT/MzPeakFile.h` is
the OpenMS interface.

## Local changes (keep when updating)

- `src/util/numpress.cpp`: includes `OpenMS/FORMAT/MSNUMPRESS/MSNumpress.h` (the copy that
  libOpenMS compiles) instead of the meson subproject's `MSNumpress.hpp`.
- `src/util/planner.cpp`: the unused `#include <print>` is dropped (libstdc++ 13 has no
  `<print>`).
- Writer: `SpectrumData` / `Util::SpectrumMetaRow` carry a scan-level ion mobility
  (`ion_mobility`, `ion_mobility_type`), written to the scans facet as the
  `ion_mobility_value` / `ion_mobility_type` columns that the reader already joins
  (`include/mzpeak/writer.h`, `include/mzpeak/util/parquet_writer.h`,
  `src/util/parquet_writer.cpp`, `src/writer.cpp`).
- Writer: `Util::minutes_from_seconds()` picks the stored minutes so that the reader's
  `* 60.0` returns the seconds exactly (`src/util/parquet_writer.cpp`, used for spectra,
  chromatograms and wavelength spectra).
- Writer: the points of a spectrum (and of a chromatogram / wavelength spectrum) are sorted along
  their axis with a stable sort, so that peaks of equal m/z keep the caller's order (`src/writer.cpp`,
  three places); the unstable sort permuted them and changed a consumer's deisotoping.
- Writer: a run whose spectra all have no points gets an empty peaks table, so that the reader
  reports the spectra (`RunArchiveWriter::finish()` in `src/writer.cpp`).
- Reader: a spectrum (or wavelength spectrum) whose metadata row has no `time` takes the
  retention time of its scan (`src/util/metadata_model.cpp`: `value_or(from_scan)` made the
  comparison a no-op, so a scan-only retention time was never assigned; found by the review).
- Reader: `Spectrum::ion_mobility_array_type()` tells the array type of the per-peak ion mobility
  array (`include/mzpeak/spectrum.h`, `src/spectrum.cpp`); a timed scan replaces an untimed
  representative scan (`src/util/metadata_model.cpp`).
- `include/mzpeak/util/compat.h`, `include/mzpeak/util/query.h`, `src/util/numpress.cpp`: `MzPeak::Util::function_ref`
  (in compat.h) replaces `boost::compat::function_ref` (Boost 1.86) in the query evaluation callbacks and
  `std::move_only_function` in the numpress decoders, whose fallback for libc++ without it was
  `boost::compat::move_only_function` (Boost 1.89). The library then needs no newer Boost than OpenMS (1.81) and no
  Boost.Compat.
- `src/writer.cpp`: `ZIP_LENGTH_TO_END` is defined as -1 where libzip does not define it (before 1.10.1, e.g.
  1.7.3 of Ubuntu 24.04), which reads the file to its end there as well.
- `src/boost_json_src.cpp` (new): Boost.JSON compiled header-only, so libOpenMS links no
  compiled Boost library.
- `sources.cmake`, `CMakeLists.txt`, this file.

## Updating

Copy `include/` and `src/` from the library, re-apply the changes above, regenerate
`sources.cmake` from the `lib_sources` list of its `meson.build`, and record the commit here.
