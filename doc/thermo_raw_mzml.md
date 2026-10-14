# Thermo RAW to mzML with metadata preservation

`ThermoRawFile` requires openms-thermo-bridge 1.0.0, the revision pinned in
`cmake/cmake_findExternalLibs.cmake` (native and managed components from the same
revision). It reads every scan with the bridge's typed single-pass API: one call per
scan returns the header (with the vendor-formatted trailer values) and, for scans that are
converted, the m/z-sorted peaks with their charges and noise arrays; no per-scan JSON is
parsed. With more than one OpenMP thread (TOPP `-threads`), up to three worker threads
read the following scans ahead on their own bridge accessors while the calling thread
converts them; the output does not depend on the thread count. It loads source SHA-1, the complete creation timestamp,
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
set `centroid = true` on both conversion paths. `FileHandler::loadExperiment()`,
through which most TOPP tools read `.raw` input, sets `centroid = true`, as
FileConverter does unless `-RawToMzML:no_peak_picking` is given. Optional
charge/noise/detector exports are disabled by default. Trailer preservation,
methods and SHA-1 are on in `Options`; each can be disabled separately to control
output size or hashing cost. The SHA-1 reads the whole file once more (5-15% of the
load time), so FileConverter computes it only with `-raw_checksum`; otherwise the
mzML `sourceFile` carries an empty SHA-1.

## Loading part of a file

`Options` also selects which scans are read. Scans that are not selected are
not decoded:

| Option | Selection |
| --- | --- |
| `first_scan`, `last_scan` | Scan number range (inclusive; 0 = unbounded) |
| `rt_min`, `rt_max` | Retention time range in seconds, `rt_min <= RT < rt_max` |
| `ms_levels` | MS levels to load (empty = all; PDA spectra have level 0) |
| `precursor_mz_min`, `precursor_mz_max` | MSn scans whose first precursor m/z is in `[min, max)`; MS1 is kept |
| `metadata_only` | Spectra without peaks, charge or noise arrays |
| `skip_chromatograms` | No TIC or detector chromatograms |

A selected MSn scan is the same as in a full load: its precursors, including the
precursor intensity that is estimated from the parent scan, do not depend on
whether the parent scan is selected (the reader decodes the parent internally if
needed). With `metadata_only`, values derived from the peaks are absent: the
precursor intensity, `lowest observed m/z` / `highest observed m/z`, and the base
peak of centroided data (the vendor scan statistics base peak is reported
instead). The TIC chromatogram always covers the whole run.

`FileHandler::loadExperiment()` passes the `PeakFileOptions` MS levels, RT and
precursor m/z ranges, metadata-only, fill-data (off = no peaks) and
skip-chromatograms settings on to the reader, so for example a search engine that
loads MS2 spectra only never decodes the MS1 scans. The m/z and intensity ranges
are applied after loading. `FileHandler` does not compute the SHA-1 of the RAW
file. Instrument configuration ids (`thermo_ic_N`) are assigned in acquisition
order over all scans read, so they are the same whichever scans are selected.

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
peaks does not reorder, truncate or reduce precision in these arrays. Scans
without available vendor noise measurements carry none of the three arrays;
absence means unavailable, not zero noise. In particular, the segmented mass
coordinates returned for scans without a native centroid stream are not exported
as noise coordinates when their corresponding noise and baseline values are absent.

## Preserving noise levels for reporter-ion quantification

Noise export is already supported by the in-process reader and mzML reader/writer.
To retain the vendor noise levels when converting a RAW file, enable it explicitly:

```sh
FileConverter -in input.raw -out output.mzML \
  -RawToMzML:reader inprocess -RawToMzML:include_noise
```

The output contains `MS:1002743` (sampled noise m/z array), `MS:1002744`
(sampled noise intensity array), and `MS:1002745` (sampled noise baseline array)
when the vendor provides them. `MzMLFile::load()` restores each array as a
double-list spectrum meta value with the corresponding name from the table above;
storing the experiment again retains these CV-annotated binary arrays.
They are not per-peak S/N values. Associating the noise grid with reporter ions
and using S/N instead of intensity for TMT quantification is a separate algorithmic
feature; exporting these arrays does not change existing quantification.

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

## Streaming

`ThermoRawFile::transform(path, consumer)` gives the same spectra and
chromatograms as `load()` with the same `Options`, but passes them one at a time
to an `Interfaces::IMSDataConsumer` instead of building an `MSExperiment`. A
first pass reads the scan headers only (no peaks). The consumer then receives:

1. `setExpectedSize()` with the exact spectrum and chromatogram counts,
2. `setExperimentalSettings()` with all instrument configurations,
3. `consumeSpectrum()` for every selected MS scan in acquisition (RT) order,
4. `consumeSpectrum()` for PDA spectra (`all_detectors`), RT-sorted among
   themselves but after the MS scans (`load()` sorts them into the MS scans),
5. `consumeChromatogram()` for the TIC, then the detector chromatograms.

Memory does not grow with the file size. For the precursor intensity estimate,
only the most recently converted parent scans are kept in a small bounded cache;
a parent that is not cached is decoded again. Detector data (`all_detectors`) is
read in one go. The SHA-1 is only computed if `Options::checksum` is set.

| Entry point | Behaviour |
| --- | --- |
| `FileHandler::transformExperiment()` | Streams mzML, mzXML and .raw into a consumer and applies the `PeakFileOptions` (MS levels, RT/m/z/intensity/precursor ranges, metadata only, skip chromatograms) as `loadExperiment()` does, for all three formats. For mzML/mzXML, `setExpectedSize()` reports the counts from the file before any selection; for .raw the scan-level options go to the reader, which reports the selected counts, and the m/z and intensity ranges are applied per spectrum |
| `FileConverter -process_lowmemory` | Streams .raw to mzML and cachedMzML only; for any other output type it logs a warning and loads the whole file |
| `OpenSwathWorkflow` | Streams .raw into the SWATH consumers; the run metadata is taken from a metadata-only read |
| `SageAdapter`, `CometAdapter` | Stream .raw into the temporary mzML |

A streamed mzML differs from an in-memory conversion only in the file content
term (the streaming writer takes it from the first spectrum) and the index
offsets.

## Building the paired branches locally

The bridge source revision is pinned in `cmake_findExternalLibs.cmake` and the
vcpkg overlay. The paired bridge is published in
https://github.com/OpenMS/openms-thermo-bridge/pull/13 (1.0.0). To use a local checkout,
use CMake's standard
`-DFETCHCONTENT_SOURCE_DIR_OPENMSTHERMOBRIDGE=/absolute/path/to/openms-thermo-bridge`
override to use the paired checkout. A matching prebuilt managed directory can
be supplied with `OPENMS_THERMO_BRIDGE_PREBUILT_MANAGED_DIR`; otherwise a .NET 8
SDK is required. Older system/vcpkg installations of the bridge are rejected by
the minimum-version check. The vcpkg overlay builds the matching managed component
from source and requires a .NET 8 SDK on the build host; it does not depend on a
pre-built binary asset. The overlay fetches the pinned GitHub source archive with
a verified SHA-512 checksum.
