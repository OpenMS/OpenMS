Vendor Formats
==============

OpenMS reads two vendor formats directly, without a conversion to {term}`mzML` first: Thermo
Fisher `.raw` files and Bruker timsTOF `.d` directories. The OpenMS installers and the pyOpenMS
wheels include both readers. For all other vendor formats, convert the data to mzML first, for
example with ProteoWizard's `msconvert`.

## Thermo Fisher `.raw`

OpenMS reads `.raw` files with Thermo's RawFileReader library, which runs on .NET. It needs
the **.NET 8 runtime** (or newer). How to install it, and how to point OpenMS to a runtime in a
non-standard location with `DOTNET_ROOT`, is described in *Reading Thermo Fisher RAW files* on
the installation page for [Linux](/about/installation/installation-on-gnu-linux.md#reading-thermo-fisher-raw-files),
[macOS](/about/installation/installation-on-macos.md#reading-thermo-fisher-raw-files) and
[Windows](/about/installation/installation-on-windows.md#reading-thermo-fisher-raw-files).

RawFileReader reading tool. Copyright © 2016 by Thermo Fisher Scientific, Inc. All rights
reserved.

### Converting with FileConverter

[FileConverter](https://archive.openms.de/openms/Documentation/nightly/latest/html/TOPP_FileConverter.html)
converts `.raw` files with one of two readers, chosen with `-RawToMzML:reader`:

| | `inprocess` | `external` |
| --- | --- | --- |
| Default on | Linux and macOS | Windows |
| Reader | built into OpenMS | ThermoRawFileParser, run as a separate program |
| Needs | .NET 8 | on Windows the .NET Framework that Windows includes; on Linux and macOS mono, which the OpenMS packages do not include |
| Output | every format FileConverter writes | mzML only |

The Windows installer puts ThermoRawFileParser on the `PATH`, so the external reader works
there without further setup. Elsewhere, give its location with `-RawToMzML:ThermoRaw_executable`
and the mono executable with `-RawToMzML:NET_executable`.

Both readers apply Thermo's peak picking by default, so the output holds centroided spectra.
`-RawToMzML:no_peak_picking` keeps the spectra as they were acquired, and
`-RawToMzML:include_noise` adds the noise data.

```bash
FileConverter -in sample.raw -out sample.mzML
```

### Reading `.raw` in other tools

These tools also accept `.raw` files as input:
CometAdapter, FeatureFinderCentroided, FeatureFinderIdentification, FeatureFinderLFQ,
FeatureFinderMetabo, FeatureFinderMetaboIdent, FeatureFinderMultiplex, FLASHDeconv,
IsobaricAnalyzer, MassTraceExtractor, MetaboliteSpectralMatcher, MS1LabeledWorkflow,
NucleicAcidSearchEngine, OpenPepXL, OpenSwathPeakMapExtractor, OpenSwathWorkflow, ProSE,
ProteomicsLFQ, SageAdapter and SimpleSearchEngine.

They read with the built-in reader, on every platform, so they need the .NET 8 runtime on
Windows as well. Like FileConverter, they apply Thermo's peak picking, so they process
centroided spectra. To process the spectra as they were acquired, which for most Orbitrap
methods means profile MS1 spectra, convert the `.raw` file with
`FileConverter -RawToMzML:no_peak_picking` first and use the mzML.

```{note}
OpenNuXL also accepts `.raw` files. It converts them with ThermoRawFileParser, like the
external reader of FileConverter, so on Linux and macOS it needs mono.
```

## Bruker timsTOF `.d`

OpenMS reads the `.d` directories of timsTOF instruments (TDF format, with an `analysis.tdf`
file) with the [opentims](https://github.com/michalsta/opentims) library. A `.d` directory
packed into a ZIP archive (`.d.zip`) can be given instead; it is unpacked into a temporary
directory while the tool runs, so that directory needs space for the unpacked data.

Bruker's SDK is not needed: OpenMS computes m/z and ion mobility (1/K0) from the calibration
stored in the file. If Bruker's SDK library (`timsdata.dll` or `libtimsdata.so`) is
installed, set the environment variable `OPENMS_BRUKER_SDK_PATH` to the library's path to
convert ion mobility with Bruker's own code instead.

These tools accept `.d` and `.d.zip` input:
CometAdapter, FeatureFinderIdentification, FeatureFinderLFQ, FeatureFinderMetaboIdent,
FileConverter, IonMobilityBinning, MetaboliteSpectralMatcher, NucleicAcidSearchEngine,
OpenSwathPeakMapExtractor, OpenSwathWorkflow, PeakPickerIM, ProSE, ProteomicsLFQ,
SageAdapter, SimpleSearchEngine and TransitionListEvidenceFilter.

The acquisition type is detected from the file:

- **DDA-PASEF**: one MS1 spectrum per frame, which stores the ion mobility of every peak, and
  one MS2 spectrum per precursor with a single ion mobility value.
- **DIA-PASEF**: MS1 as above, and one MS2 spectrum per frame and isolation window, which
  stores the ion mobility of every peak.

FileConverter, CometAdapter and PeakPickerIM have `bruker:` options to change this, for
example `-bruker:export_mode` to force the per-precursor layout (`spectrum`) or raw frames
(`frame`). FileConverter can also skip the MS1 spectra (`-bruker:load_ms1 false`),
recalibrate m/z, aggregate neighbouring frames and centroid along the ion mobility axis; see
its documentation. ProteomicsLFQ
centroids the MS1 frames along the ion mobility axis before feature detection. The other
tools read with the default settings.

To convert a `.d` directory to mzML:

```bash
FileConverter -in sample.d -out sample.mzML
```

## pyOpenMS

pyOpenMS reads both formats as well; see
[Vendor formats](https://pyopenms.readthedocs.io/en/latest/user_guide/vendor_formats.html) in
the pyOpenMS user guide.
