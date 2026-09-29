File Formats
============

OpenMS reads spectra from {term}`mzML` and stores its results in its own XML formats: featureXML for features,
consensusXML for features linked across runs, and idXML for identifications. This page describes three additions of
OpenMS 3.6: compressed input and output, Parquet versions of the result formats, and imaging data. For Thermo `.raw`
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
