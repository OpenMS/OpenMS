# Thermo RAW to mzML with metadata preservation

`ThermoRawFile` requires the openms-thermo-bridge revision pinned in
`cmake/cmake_findExternalLibs.cmake` (native and managed components from the same
revision). It loads source SHA-1, the complete creation timestamp,
sample fields, user sample fields, embedded method texts and the full scan trailer.
Vendor values without a confirmed standardized unit are kept under `Thermo ...`
metadata keys; they are not assigned guessed units.

```cpp
OpenMS::ThermoRawFile reader;
OpenMS::ThermoRawFile::Options options;
options.centroid = true;       // matches TRFP's default peak picking
options.charge_data = true;    // TRFP --chargeData
options.noise_data = true;     // TRFP --noiseData
options.all_detectors = true;  // TRFP --allDetectors
reader.setOptions(options);
OpenMS::MSExperiment experiment;
reader.load("input.raw", experiment);
OpenMS::MzMLFile().store("output.mzML", experiment);
```

The same options are exposed in pyOpenMS as `ThermoRawFileOptions` and through
`ThermoRawFile.getOptions()` / `setOptions()`. The default preserves the acquired
profile/centroid representation; compare with TRFP's no-peak-picking mode or
set `centroid = true` on both conversion paths. Optional charge/noise/detector
exports are disabled by default. Trailer preservation, methods and SHA-1 are on;
each can be disabled separately to control output size or hashing cost.

| Metadata | OpenMS representation |
| --- | --- |
| Default instrument | `MSExperiment.getInstrument()` |
| Other analyzer/source configurations | `getInstrumentConfigurations()`, referenced by acquisition `instrument_configuration_ref` |
| Isolation target versus monoisotopic selection | `Precursor.getMZ()` is the selected ion; `isolation window target m/z` metadata preserves the target |
| MSn hierarchy and SPS | Ordered precursor list with each precursor's `spectrum_ref`; ancestor descriptors retain their own references |
| Supplemental activation | Separate supplemental CID/HCD and collision energy metadata; main activation remains distinct |
| Full trailer | Acquisition `Thermo trailer extra`, JSON label/value pairs |
| Embedded methods | Experiment `Thermo instrument methods`, JSON string array |
| Native peak charges | Integer data array named `charge array` |
| Sampled noise | Double-list metadata named `sampled noise m/z array`, `sampled noise intensity array`, `sampled noise baseline array`; written as independent 64-bit mzML binary arrays |
| PDA spectra | Wavelength coordinates in the peak position slot, absorption scan mode and `mzml coordinate array = wavelength` / `mzml intensity array = absorption` |
| Detector signals | Chromatograms with original labels/units; pressure/flow/absorption CV array types retained |

Sampled noise is deliberately separate from peak annotations. Its point count
and mass grid can differ from the spectrum; sorting or restricting spectrum
peaks does not reorder, truncate or reduce precision in these arrays.

The target is comparable mzML content, not identical XML bytes: OpenMS writes
retention times in seconds, has its own software/data-processing records and
compression/index layout, and its peak intensity container is float32. Extra
vendor metadata is preserved beyond what TRFP writes to mzML. Match centroiding
and optional export flags before comparing the files. Independently sampled
noise values and m/z coordinates retain double precision.

Regression coverage includes the standalone `ThermoRawFileMetadata_test` for
precursor reconstruction, `MzMLFile_test` for synthetic metadata round trips,
`ExperimentalSettings_test` for configuration ownership, and `ThermoRawFile_test`
for RAW-to-mzML integration when `ENABLE_THERMO_RAW_TESTS=ON`. Python API tests
are in `test_ThermoRawFile.py`.

## Building the paired branches locally

The bridge source revision is pinned in `cmake_findExternalLibs.cmake` and the
vcpkg overlay. The paired bridge is published in
https://github.com/OpenMS/openms-thermo-bridge/pull/11. To use a local checkout,
use CMake's standard
`-DFETCHCONTENT_SOURCE_DIR_OPENMSTHERMOBRIDGE=/absolute/path/to/openms-thermo-bridge`
override to use the paired checkout. A matching prebuilt managed directory can
be supplied with `OPENMS_THERMO_BRIDGE_PREBUILT_MANAGED_DIR`; otherwise a .NET 8
SDK is required. Older system/vcpkg installations of the bridge are rejected by
the minimum-version check. The vcpkg overlay builds the matching managed component
from source and requires a .NET 8 SDK on the build host; it does not depend on a
pre-built binary asset. The overlay fetches the pinned GitHub source archive with
a verified SHA-512 checksum.
