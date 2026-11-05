File Formats
============

OpenMS reads spectra from {term}`mzML` and stores its results in its own XML formats: featureXML for features,
consensusXML for features linked across runs, and idXML for identifications. This page describes three additions of
OpenMS 3.6: compressed input and output, Parquet versions of the result formats, and imaging data. It also covers
mzPeak and reading XML files from a URL, new in OpenMS 3.7. For Thermo `.raw`
files and Bruker `.d` directories, see [Vendor formats](vendor-formats.md).

## Compressed files

### Reading

Tools read compressed XML files directly, without unpacking them first. This works for mzML, mzXML, mzData,
featureXML, consensusXML, idXML, pepXML, protXML, traML, trafoXML and INI files, among others. Add `.gz` (gzip),
`.bz2` (bzip2) or `.zip` to the file name, for example `sample.mzML.gz`. A ZIP archive must contain exactly one file.

```{note}
In OpenMS 3.6, mzIdentML files (`.mzid`) cannot be read compressed. Unpack them first.
```

Other formats, for example MGF or text files, have to be unpacked before a tool can read them; a tool refuses a
compressed file of such a format. A Bruker `.d` directory packed into a ZIP archive (`.d.zip`) is described in
[Vendor formats](vendor-formats.md).

### Writing

mzML, mzXML, mzData, featureXML, consensusXML, traML and mzIdentML output is compressed when the output file name
ends in `.gz` (gzip) or `.bz2` (bzip2):

```bash
FileConverter -in sample.mzML -out sample.mzML.gz
```

The other formats, idXML and trafoXML among them, cannot be written compressed, and ZIP output is not supported.
A tool refuses such an output name with an error rather than write an uncompressed file under it: give the output
an uncompressed name, and compress the file afterwards if needed. The same applies to the modes that write spectra
one by one, for example FileConverter with `-process_lowmemory`: they write uncompressed mzML only.

### Compressing the peak data inside mzML

The peak data inside an mzML file can be compressed as well. With `-zstd_compression`, new in OpenMS 3.6,
FileConverter compresses it with Zstandard (zstd); FileFilter has the same option as
`-peak_options:zstd_compression true`:

```bash
FileConverter -in sample.mzML -out sample_zstd.mzML -zstd_compression
```

zstd is lossless and usually gives smaller files than zlib, and it can be combined with the lossy numpress compression
of `-lossy_compression`. Tools read such files like any other mzML file. Not every program outside OpenMS reads
zstd-compressed mzML yet, so leave the option off for files that other software has to read.

## Parquet bundles

OpenMS 3.6 can store identifications, feature maps and consensus maps as [Apache Parquet](https://parquet.apache.org/)
tables instead of idXML, featureXML and consensusXML. Such a *bundle* is a directory of Parquet files:

| Format | Instead of | Files in the directory |
| --- | --- | --- |
| `.idparquet` | idXML | `psms.parquet`, `proteins.parquet`, `protein_groups.parquet`, `search_params.parquet` |
| `.featureparquet` | featureXML | `features.parquet` and the four files of `.idparquet` |
| `.consensusparquet` | consensusXML | `consensus_features.parquet` and the four files of `.idparquet` |

Any Parquet reader can load these tables without OpenMS, for example pandas or pyarrow in Python, the arrow package
in R, or DuckDB.

These tools read or write the bundles:

- `.idparquet`: CometAdapter, FeatureFinderIdentification, IDFileConverter, IDFilter, IDMerger, IDRipper,
  IDScoreSwitcher, IsobaricWorkflow, MapAlignerIdentification, MapRTTransformer, MS1LabeledWorkflow,
  MSGFPlusAdapter, MzTabExporter, PeptideIndexer, PercolatorAdapter, PSMFeatureExtractor, SageAdapter and
  TextExporter.
- `.featureparquet`: FeatureFinderIdentification, the FeatureLinker tools, FileConverter, IDConflictResolver,
  MapAlignerIdentification, MapRTTransformer, MzTabExporter and TextExporter.
- `.consensusparquet`: the FeatureLinker tools, FileConverter, IDConflictResolver, IDFilter,
  MapAlignerIdentification, MapRTTransformer, MzTabExporter and TextExporter.

To convert between the XML formats and the bundles, use FileConverter or ParquetConverter for feature and consensus
maps, and IDFileConverter for identifications:

```bash
FileConverter -in sample.featureXML -out sample.featureparquet
IDFileConverter -in sample.idXML -out sample.idparquet
```

ParquetDiff compares two Parquet files, for example the same table of two bundles. It matches rows by a key and
compares values with a tolerance, so a different row order or tiny floating-point differences do not count as
changes.

## Imaging data (imzML)

OpenMS 3.6 reads and writes imzML, the standard format for mass spectrometry imaging, in its library and in pyOpenMS.
The TOPP tools and TOPPView do not accept imzML files yet.

## mzPeak

OpenMS 3.7 reads and writes [mzPeak](https://github.com/HUPO-PSI/mzPeak) (HUPO-PSI, version 0.9.0), a spectrum format
that stores the data as Parquet tables. An mzPeak run is a ZIP archive (`.mzpeak`) of Parquet tables, or a directory
with the same members. It is a raw-data format like mzML: FileConverter converts to and from it, and FileHandler,
FileInfo and the tools that read spectra through FileHandler recognise it:

```bash
FileConverter -in sample.mzML -out sample.mzpeak
FileConverter -in sample.mzpeak -out sample.mzML
```

mzPeak support is a build option, `OPENMS_ENABLE_MZPEAK`. It is left out of a build whose compiler and Boost lack
what the vendored mzPeak library needs, and the format is then not available. In pyOpenMS the format is
`FileType.MZPEAK`; the library class is `MzPeakFile`.

## Reading XML from a URL

OpenMS 3.7 reads XML files, mzML for example, from `http://`, `https://` and `ftp://` URLs with libcurl, instead of
from a local file. The data is streamed while it is parsed, and gzip and bzip2 content is decompressed on the fly.
Redirects are followed, but only to HTTP(S) and FTP. An HTTP error such as 404 is reported as a missing file, and a
transfer slower than 1 byte per second for two minutes fails. Credentials in a URL are removed from log and error
messages.
