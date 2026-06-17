File Handling
=============

## General information about peak and feature maps

For general information about a peak or feature map, use the `FileInfo` tool.

- It can print RT, m/z and intensity ranges, the overall number of peaks, and the distribution of MS levels.
- It can print a statistical summary of intensities.
- It can print some meta information.
- It can validate XML files against their schema.
- It can check for corrupt data in peak files See the `FileInfo –help` for details.

## Problems with input files

If you are experiencing problems while processing an XML file, check if the file does validate against the XML schema:

`FileInfo -v -in infile.mzML`

Validation is available for several file formats including mzML, mzData, mzXML, featureXML and idXML.

Another frequently-occurring problem is corrupt data. You can check for corrupt data in peak files with `FileInfo` as well:

`FileInfo -c -in infile.mzML`

## Converting your files to mzML

The TOPP tools work only on the HUPO-PSI `mzML` format. If you need to convert *mzData*, *mzXML* or *ANDI/MS* data to
*mzML*, use the FileConverter, e.g.

`FileConverter -in infile.mzXML -out outfile.mzML`

For format names as file extension, the tool derives the format from the extension. For other extensions, the file
formats of the input and output file can be given explicitly.

## Converting between DTA and mzML

Sequest DTA files can be extracted from a mzML file using the `DTAExtractor`:

`DTAExtractor -in infile.mzML -out outfile`

The retention time of a scan, the precursor mass-to-charge ratio (for MS/MS scans) and the file extension are appended
to the output file name.

To combine several files (e.g. DTA files) to an `mzML` file use the `FileMerger`:

`FileMerger -in infile_list.txt -out outfile.mzML`

The retention times of the scans can be generated, taken from the `infile_list.txt` or can be extracted from the DTA
file names. See the FileMerger documentation for details.

## Extracting part of the data from a file

To extract part of the data from an `mzML` file, use the `FileFilter` tool. It allows filtering for RT, `m/z` and
intensity range or for MS level. To extract the MS/MS scans between retention time `100` and `1500`, use the following
command:

`FileFilter -in infile.mzML -levels 2 -rt 100:1500 -out outfile.mzML`

## Conversion Between OpenMS XML Formats and Text Formats

### Export of OpenMS XML formats

As TOPP offers no functionality for statistical analysis, this step is normally done using external statistics packages.
In order to export the OpenMS XML formats into an appropriate format for these packages the TOPP **TextExporter** can be
used.

It converts the the following OpenMS XML formats to text files:

- featureXML
- idXML
- consensusXML

The use of the `TextExporter` is is very simple:

`TextExporter -in infile.idXML -out outfile.txt`

### Import of feature data to OpenMS

OpenMS offers a lot of visualization and analysis functionality for feature data.
Feature data in text format, e.g. from other analysis tools, can be imported using the `TextImporter`. The default
mode accepts comma separated values containing the following columns: RT, m/z, intensity. Additionally meta data
columns may follow. If meta data is used, meta data column names have to be specified in a header line. Without headers:

```bash
1201	503.123	1435000
1201	1006.246	1235200
```

Or with headers:

```bash
RT	m/z	Int	isHeavy	myMeta
1201	503.123	1435000	true	2
1201	1006.246	1235200	maybe	1
```

Example invocation:

`TextImporter -in infile.txt -out outfile.featureXML`

The tool also supports data from msInspect,SpecArray and Kroenik(Hardkloer sibling), just specify the `-mode` option
accordingly.

### Import of protein/peptide identification data to OpenMS

Peptide/protein identification data from several identification engines can be converted to idXML format using the
**IDFileConverter** tool.

It can currently read the following formats:
- Sequest output folder
- pepXML file
- idXML file

It can currently write the following formats:

- pepXML
- idXML

This example shows how to convert pepXML to idXML:

`IDFileConverter -in infile.pepXML -out outfile.idXML`