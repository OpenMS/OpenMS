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

## Random access

`OnDiscThermoRawExperiment` keeps the RAW file open and reads single spectra on
request, like `OnDiscMSExperiment` does for indexed mzML. Both implement
`IMSExperimentRandomAccess`. `openFile(path, options)` runs the header pass of
`transform()` once. It builds an index of the selected scans (scan number,
RT, MS level, native id) and resolves the precursors in acquisition order. After
that, `getSpectrum(i)` converts any scan on its own and returns the `i`-th spectrum
of `load()` with the same `Options`. This includes the precursor intensity
estimate, which decodes the parent MS1 scan on demand and keeps recent parents in
a small LRU cache. The `Options` (MS levels, scan/RT/m/z ranges, peak picking,
`metadata_only`, `skip_chromatograms`, `all_detectors`) define which spectra are
indexed. `Options::checksum` is ignored.

| Method | Returns |
| --- | --- |
| `getNrSpectra()`, `getNrChromatograms()` | Counts as in `load()` |
| `getSpectrum(i)`, `getSpectrumMetadata(i)` | Spectrum with peaks, or without peaks (no decoding) |
| `getSpectrumRT(i)`, `getSpectrumMSLevel(i)`, `getSpectrumNativeID(i)` | Index lookups (no file access) |
| `getScanNumber(i)`, `getSpectrumIndexByScanNumber(scan)`, `getSpectrumByNativeId(id)` | Scan number mapping |
| `findFirstRT(rt)` | Index of the first spectrum with RT >= `rt` |
| `getChromatogram(i)` | TIC and detector chromatograms |
| `getExperimentalSettings()` | Settings as in `load()` |

Reads on one object are serialized. For parallel reads, give each thread its own
`clone()`. Clones share the immutable index and each get a bridge accessor of their
own, so they read concurrently and return the same spectra as sequential reads.
Concurrent reads need an OpenMP build of OpenMS, because the meta value registry is
guarded by OpenMP critical sections. PDA spectra (`all_detectors`) are read
completely at open and kept in memory.

`SpectrumAccessOnDisc` exposes any `IMSExperimentRandomAccess` as an
`OpenSwath::ISpectrumAccess`. It can show all spectra or an RT-sorted subset of
indices, for example one SWATH window. `lightClone()` clones the experiment, so
clones can be used from different threads.

`OpenSwathWorkflow -readOptions onDisc` uses these classes for .raw input
(and `OnDiscMSExperiment` for indexed mzML; `SwathFile::loadOnDisc()` chooses
the backend by file type). The
spectra are sorted into SWATH windows from the index, and each window reads its
spectra from the file during extraction. Nothing is loaded or cached in advance.
The results are the same as with `normal`. Two restrictions apply:

- OpenSwathWorkflow does not process FAIMS compensation voltages separately with
  any read option. With `onDisc`, input with more than one CV is rejected instead
  of merging the CVs into one set of SWATH windows; single-CV input works.
  `SwathFile::loadOnDiscByFAIMSCV()` splits MS1
  and MS2 spectra by CV, as `IMDataConverter::splitByFAIMSCV()` does, and
  builds the SWATH maps within each CV.
- The spectra are read from the file again during extraction, so a plugin
  consumer (such as the QC of `-out_qc`) may only observe them. A consumer that
  modifies spectra raises an error; use `normal` or `cache` for it.

In pyOpenMS:

```python
import pyopenms as oms

exp = oms.OnDiscThermoRawExperiment()
exp.openFile("run.raw")          # optional second argument: ThermoRawFile.Options
spec = exp[0]                    # = exp.getSpectrum(0); negative indices work
meta = exp.getSpectrumMetadata(10)
first = exp.findFirstRT(600.0)
for spec in exp:                 # iterate in index order
    pass

class Counter:
    def __init__(self): self.n = 0
    def setExpectedSize(self, s, c): pass
    def setExperimentalSettings(self, s): pass
    def consumeSpectrum(self, s): self.n += 1
    def consumeChromatogram(self, c): pass

oms.ThermoRawFile().transform("run.raw", Counter())
```

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
