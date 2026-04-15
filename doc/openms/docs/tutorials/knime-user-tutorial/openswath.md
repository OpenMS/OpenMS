OpenSWATH
=========

## Introduction

[OpenSWATH](http://openswath.org/en/latest/index.html) [^3] allows the analysis of LC-MS/MS DIA (data independent acquisition) data using the approach described by Gillet *et al*. [^4]. The DIA approach described there uses 32 cycles to iterate through precursor ion windows from 400-426 Da to 1175-1201 Da and at each step acquires a complete, multiplexed fragment ion spectrum of all precursors present in that window. After 32 fragmentations (or 3.2 seconds), the cycle is restarted and the first window (400-426 Da) is fragmented again, thus delivering complete “snapshots” of all fragments of a specific window every 3.2 seconds.
The analysis approach described by Gillet et al. extracts ion traces of specific fragment ions from all MS2 spectra that have the same precursor isolation window, thus generating data that is very similar to SRM traces.

## Installation of OpenSWATH

OpenSWATH has been fully integrated since OpenMS 1.10 [^2]<sup>,</sup> [^1]<sup>,</sup> [^5]<sup>,</sup> [^6]<sup>,</sup> [^7].

## Installation of mProphet

mProphet[^8] is available as standalone script in {path}`External_Tools,mProphet` or can be downloaded [here](https://github.com/OpenMS/OpenMS-Tutorials/releases/download/data-and-tools-OpenMSv2.0.0/External_Tools.zip). [R](http://www.r-project.org/) and the package [MASS](http://cran.r-project.org/web/packages/MASS/) are further required to execute mProphet. Please obtain a version for either Windows, Mac or Linux directly from CRAN.
PyProphet, a much faster reimplementation of the mProphet algorithm is available from [PyPI](https://pypi.python.org/pypi/pyprophet/). The usage of pyprophet instead of mProphet is suggested for large-scale applications.

mProphet will be used in this tutorial.

## Generating the Assay Library

### Generating TraML from transition lists

OpenSWATH requires an assay library to be supplied in the TraML format[^9]. To enable manual editing of transition lists, the TOPP tool **TargetedFileConverter** is available, which uses tab separated files as input. Example datasets are provided in {path}`ExampleData,OpenSWATH,assay` or downloaded from [here](https://github.com/OpenMS/OpenMS-Tutorials/releases/download/data-and-tools-OpenMSv2.0.0/Example_Data.zip). Please note that the transition lists need to be named `.tsv`.

The header of the transition list contains the following variables (with example values in brackets):

**Required Columns**:
`PrecursorMz`

The mass-to-charge (m/z) of the precursor ion. (924.539)

`ProductMz`

The mass-to-charge (m/z) of the product or fragment ion. (728.99)

`LibraryIntensity`

The relative intensity of the transition. (0.74)

`NormalizedRetentionTime`

The normalized retention time (or iRT)[^10] of the peptide. (26.5)

**Targeted Proteomics Columns**
`ProteinId`

A unique identifier for the protein. (AQUA4SWATH_HMLangeA)

`PeptideSequence`

The unmodified peptide sequence. (ADSTGTLVITDPTR)

`ModifiedPeptideSequence`

The peptide sequence with UniMod modifications. (ADSTGTLVITDPTR(UniMod:267))

`PrecursorCharge`

The precursor ion charge. (2)

`ProductCharge`

The product ion charge. (2)

**Grouping Columns**:
`TransitionGroupId`

A unique identifier for the transition group.
(AQUA4SWATH_HMLangeA_ADSTGTLVITDPTR(UniMod:267)/2)

`TransitionId`

A unique identifier for the transition.
(AQUA4SWATH_HMLangeA_ADSTGTLVITDPTR(UniMod:267)/2_y8)

`Decoy`

A binary value whether the transition is target or decoy. (target: 0, decoy: 1)

`PeptideGroupLabel`

Which label group the peptide belongs to.

`DetectingTransition`

Use transition for peak group detection. (1)

`IdentifyingTransition`

Use transition for peptidoform inference using IPF. (0)

`QuantifyingTransition`

Use transition to quantify peak group. (1)

For further instructions about generic transition list and assay library generation please see the following [link](http://openswath.org/en/latest/docs/generic.html).
To convert transitions lists to TraML, use the TargetedFileConverter: Please use the absolute path to your OpenMS installation.

**Linux or Mac**

On the Terminal:

```bash
 TargetedFileConverter -in OpenSWATH_SGS_AssayLibrary_woDecoy.tsv -out OpenSWATH_SGS_AssayLibrary_woDecoy.TraML
```

**Windows**

On the TOPP command:

```bash
 TargetedFileConverter.exe -in OpenSWATH_SGS_AssayLibrary_woDecoy.tsv -out OpenSWATH_SGS_AssayLibrary_woDecoy.TraML
```

### Appending decoys to a TraML file

In addition to the target assays, OpenSWATH requires decoy assays in the library which are later used for classification and error rate estimation. For the decoy generation it is crucial that the decoys represent the targets in a realistic but unnatural manner without interfering with the targets. The methods for decoy generation implemented in OpenSWATH include ’shuffle’, ’pseudo-reverse’, ’reverse’ and ’shift’. To append decoys to a TraML, the TOPP tool **OpenSwathDecoyGenerator** can be used: Please use the absolute path to your OpenMS installation.

**Linux or Mac**

On the Terminal:

```bash
OpenSwathDecoyGenerator -in OpenSWATH_SGS_AssayLibrary_woDecoy.TraML -out OpenSWATH_SGS_AssayLibrary.TraML -method shuffle -switchKR false
```

**Windows**

On the TOPP command:

```bash
OpenSwathDecoyGenerator.exe -in OpenSWATH_SGS_AssayLibrary_woDecoy.TraML -out OpenSWATH_SGS_AssayLibrary.TraML -method shuffle -switchKR false
```

## OpenSWATH KNIME

An example KNIME workflow for OpenSWATH is supplied in `Workflows` (<a href="#figure-44">Fig. 44</a>). The example dataset can be used for this workflow (filenames in brackets):

1. Open {path}`Workflows,OpenSWATH.knwf` in KNIME: **File** > **Import KNIME Workflow...**
2. Select the normalized retention time (iRT) assay library in TraML format by double-clicking on node `File Importer` > **iRT Assay Library**. ({path}`ExampleData,OpenSWATH,assay,OpenSWATHiRTAssayLibrary.TraML`).
3. Select the SWATH MS data in mzML format as input by double-clicking on node **Input File** > **SWATH-MS files**. ({path}`ExampleData,OpenSWATH,data,splitnapedroL120420x010SW-*.nf.pp.mzML`).
4. Select the target peptide assay library in TraML format as input by double-clicking on node `Input Files` > **Assay Library**. ({path}`ExampleData,OpenSWATH,assay,OpenSWATHSGSAssayLibrary.TraML`).
5. Set the output destination by double-clicking on node `Output File`.
6. Run the workflow.

The resulting output can be found at your selected path, which will be used as input for mProphet. Execute the script on the Terminal (Linux or Mac) or cmd.exe (Windows) in {path}`ExampleData,OpenSWATH,result`. Please use the absolute path to your R installation and the result file:

```r
R --slave --args bin_dir=../../../External_Tools/mProphet/ mquest=OpenSWATH_quant.tsv workflow=LABEL_FREE num_xval=5 run_log=FALSE write_classifier=1 write_all_pg=1 < ../../../External_Tools/mProphet/mProphet.R
```

or for Windows:

```r
"C:\Program Files\R\R-3.5.1\bin\x86\R.exe" --slave --args bin_dir=../../../External_Tools/mProphet/ mquest=OpenSWATH_quant.tsv workflow=LABEL_FREE num_xval=5 run_log=FALSE write_classifier=1 write_all_pg=1 < ../../../External_Tools/mProphet/mProphet.R
```

The main output will be called: {path}`OpenSWATH,result,mProphetxallxpeakgroups.xls` with statistical information available in {path}`OpenSWATH,result,mProphet.pdf`.

Please note that due to the semi-supervised machine learning approach of mProphet the results differ slightly when mProphet is executed several times.

(Figure_44)=
|![OpenSWATH KNIME Workflow.](/_images/openms-user-tutorial/openswath/OpenSWATHWF.png)|
|:--:|
|Figure 44: OpenSWATH KNIME Workflow.|

Additionally, the chromatogram output (.mzML) can be visualized for inspection with TOPPView.
For additional instructions on how to use pyProphet instead of mProphet please have a look at the [PyProphet Legacy Workflow](http://openswath.org/en/latest/docs/pyprophet_legacy.html). If you want to use the SQLite-based workflow in your lab in the future, please have a look [here](http://openswath.org/en/latest/docs/pyprophet.html). The SQLite-based workflow will not be part of the tutorial.

## From the example dataset to real-life applications

The sample dataset used in this tutorial is part of the larger SWATH MS Gold Standard (SGS) dataset which is described in the publication of Roest *et al.*[^3]. It contains one of 90 SWATH-MS runs with significant data reduction (peak picking of the raw, profile data) to make file transfer and working with it easier. Usually SWATH-MS datasets are huge with several gigabyte per run. Especially when complex samples in combination with large assay libraries are analyzed, the TOPP tool based workflow requires a lot of computational resources. Additional information and instruction can be found at the following [link](http://openswath.org/en/latest/).


## References

[^1]: M. Sturm, A. Bertsch, C. Gröpl, A. Hildebrandt, R. Hussong, E. Lange, N. Pfeifer,
O. Schulz-Trieglaff, A. Zerck, K. Reinert, and O. Kohlbacher, <a href="http://dx.doi.org/10.1186/1471-2105-9-163">OpenMS - an opensource software framework for mass spectrometry</a>., BMC bioinformatics 9(1)
(2008), <a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-163">doi:10.1186/1471-2105-9-163</a>. 7, 83

[^2]: O. Kohlbacher, K. Reinert, C. Gröpl, E. Lange, N. Pfeifer, O. Schulz-Trieglaff, and
M. Sturm, <a href="http://view.ncbi.nlm.nih.gov/pubmed/17237091">TOPP–the OpenMS proteomics pipeline</a>., Bioinformatics 23(2) (Jan.
2007). 7, 83

[^3]: H. L. Röst, G. Rosenberger, P. Navarro, L. Gillet, S. M. Miladinovic, O. T. Schubert, W. Wolski, B. C. Collins, J. Malmstrom, L. Malmström, and R. Aebersold,
OpenSWATH enables automated, targeted analysis of data-independent acquisition MS data, Nature Biotechnology 32(3), 219–223 (Mar. 2014). 83, 87

[^4]: L. C. Gillet, P. Navarro, S. Tate, H. Röst, N. Selevsek, L. Reiter, R. Bonner, and
R. Aebersold, Targeted Data Extraction of the MS/MS Spectra Generated by Data-independent Acquisition: A New Concept for Consistent and Accurate Proteome Analysis, Molecular & Cellular Proteomics 11(6) (June 2012), doi:10.1074/mcp.O111.016717, 83

[^5]: A. Bertsch, C. Gröpl, K. Reinert, and O. Kohlbacher, <a href="https://link.springer.com/protocol/10.1007/978-1-60761-987-1_23">OpenMS and TOPP: open source software for LC-MS data analysis.</a>, Methods in molecular biology (Clifton,
N.J.) 696, 353–367 (2011), <a href="https://link.springer.com/protocol/10.1007/978-1-60761-987-1_23">doi:10.1007/978-1-60761-987-1_23</a>. 83

[^6]: H. L. Röst, T. Sachsenberg, S. Aiche, C. Bielow, H. Weisser, F. Aicheler, S. Andreotti, H.-c. Ehrlich, P. Gutenbrunner, E. Kenar, X. Liang, S. Nahnsen, L. Nilse,
J. Pfeuffer, G. Rosenberger, M. Rurik, U. Schmitt, J. Veit, M. Walzer, D. Wojnar,
W. E. Wolski, O. Schilling, J. S. Choudhary, L. Malmström, R. Aebersold, K. Reinert, and O. Kohlbacher, <a href="https://www.nature.com/articles/nmeth.3959">OpenMS: a flexible open-source software platform
for mass spectrometry data analysis</a>, Nat. Methods 13(9), 741–748 (sep 2016),<a href="https://www.nature.com/articles/nmeth.3959">doi:10.1038/nmeth.3959</a>. 83

[^7]: J. Pfeuffer, T. Sachsenberg, O. Alka, M. Walzer, A. Fillbrunn, L. Nilse, O. Schilling,
K. Reinert, and O. Kohlbacher, OpenMS - A platform for reproducible analysis
of mass spectrometry data, J. Biotechnol. 261(February), 142–148 (2017), <a href="https://www.sciencedirect.com/science/article/pii/S0168165617302511?via%3Dihub">doi:10.1016/j.jbiotec.2017.05.016</a>. 83

[^8]: L. Reiter, O. Rinner, P. Picotti, R. Huttenhain, M. Beck, M.-Y. Brusniak, M. O. Hengartner, and R. Aebersold, <a href="https://www.nature.com/articles/nmeth.1584">mProphet: automated data processing and statistical validation for large-scale SRM experiments</a>, Nature Methods 8(5), 430–435 (May
2011), <a href="https://www.nature.com/articles/nmeth.1584">doi:10.1038/nmeth.1584</a>. 83

[^9]: E. W. Deutsch, M. Chambers, S. Neumann, F. Levander, P.-A. Binz, J. Shofstahl,
D. S. Campbell, L. Mendoza, D. Ovelleiro, K. Helsens, L. Martens, R. Aebersold,
R. L. Moritz, and M.-Y. Brusniak, TraML—A Standard Format for Exchange of Selected Reaction Monitoring Transition Lists, Molecular & Cellular Proteomics
11(4) (Apr. 2012), doi:10.1074/mcp.R111.015040, 84

[^10]: C. Escher, L. Reiter, B. MacLean, R. Ossola, F. Herzog, J. Chilton, M. J. MacCoss, and O. Rinner, <a href="https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/pmic.201100463">Using iRT, a normalized retention time for more targeted measurement of peptides</a>., Proteomics 12(8), 1111–1121 (Apr. 2012),
<a href="https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/pmic.201100463">doi:10.1002/pmic.201100463</a>. 84
