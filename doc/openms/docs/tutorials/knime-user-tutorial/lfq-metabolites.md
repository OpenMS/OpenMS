Label-free quantification of metabolites
========================================

## Introduction

Quantification and identification of chemical compounds are basic tasks in metabolomic studies. In this tutorial session we construct a UPLC-MS based, label-free quantification and identification workflow. Following quantification and identification we then perform statistical downstream analysis to detect quantification values that differ significantly between two conditions. This approach can, for example, be used to detect biomarkers. 

Here, we analyze a dataset derived from bacterial cytosolic fractions to investigate the metabolic effects of fosfomycin, an antibiotic that inhibits a key step in peptidoglycan biosynthesis. The study is based on *Bacillus subtilis* cultures subjected to different treatment conditions.

The dataset consists of four sample types:

- **Blank** – A sample without any analyte
- **Control** – An untreated *B. subtilis* sample, representing baseline metabolic activity.
- **Treatment** – A *B. subtilis* sample treated with fosfomycin
- **Pool** – A composite sample combining aliquots from the three other conditions

The primary goal of this analysis is to evaluate the metabolic response of *B. subtilis* to fosfomycin treatment, with a specific focus on its effects on peptidoglycan biosynthesis pathways.

## Basics of non-targeted metabolomics data analysis

For the metabolite quantification we choose an approach similar to the one used for peptides, but this time based on the OpenMS `FeatureFinderMetabo` method. This feature finder again collects peak picked data into individual mass traces. The reason
why we need a different feature finder for metabolites lies in the step after trace detection: the aggregation of isotopic traces belonging to the same compound ion into the same feature. Compared to peptides with their averagine model, small molecules
have very different isotopic distributions. To group small molecule mass traces correctly, an aggregation model tailored to small molecules is thus needed.

- Create a new workflow called for instance ”Metabolomics”.
- Add an `File Importer` node and configure it with one mzML file [from the](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Tutorials/Example_Data/Metabolomics/datasets/) {path}`Example_Data,Metabolomics,datasets`.
- Add a `FeatureFinderMetabo` node (from **Community Nodes** > **OpenMS** > **Quantitation**) and
  connect the first output port of the `File Importer` to the `FeatureFinderMetabo`.
-  For an optimal result adjust the following settings. Please note that some of these are advanced parameters.
- Connect a `Output Folder` to the output of the `FeatureFinderMetabo` (see <a href="#figure-27">Fig. 27</a>).

(Figure_27)=
|![FeatureFinderMetabo workflow](/_images/openms-user-tutorial/metabo/minimal_FFM_wf.png)|
|:--:|
|Figure 27: FeatureFinderMetabo workflow|

In the following advanced parameters will be highlighted. These parameter can be altered if the `Show advanced parameter` field in the specific tool is activated.

|**parameter**|**value**|
|:------------|:--------|
|*algorithm*→*common*→*noise_threshold_int*|1000.0|
|*algorithm*→*common*→*chrom_peak_snr*|3.0|
|*algorithm*→*common*→*chrom_fwhm*|5.0|
|*algorithm*→*mtd*→*mass_error_ppm*|10.0|
|*algorithm*→*ffm*→*isotope_filtering_model*|'none'|
|*algorithm*→*ffm*→*remove_single_traces*|true|
|*algorithm*→*ffm*→*report_convex_hulls*|true|

The parameters change the behavior of `FeatureFinderMetabo` as follows:

- **noise_threshold_int**: Intensity threshold below which peaks are regarded as noise.
- **chrom_peak_snr**: Minimum signal-to-noise a mass trace should have.
- **chrom_fwhm**: The expected chromatographic peak width in seconds.
- **mass_error_ppm**: Allowed mass deviation (in ppm)
- **isotope_filtering_model**: Remove/score candidate assemblies based on isotope intensities. SVM isotope models for metabolites were trained with either 2% or 5% RMS error. For peptides, an averagine cosine scoring is used. Select the appropriate noise model according to the quality of measurement or MS device. (valid: 'metabolites (2% RMS)', 'metabolites (5% RMS)', 'peptides', 'none')
- **remove_single_traces**: If set to true, unassembled traces are removed (single traces).
- **report_convex_hulls**: If set to true, convex hulls including mass traces will be reported for all identified features. This increases the output size considerably.

The output file .featureXML can be visualized with TOPPView on top of the used `.mzML` file - in a so called layer - to look at the identified features.

First start TOPPView and open the example `.mzML` file (see <a href="#figure-28">Fig. 28</a>). Afterwards open the `.featureXML` output as new layer (see <a href="#figure-29">Fig. 29</a>). The overlay is depicted in <a href="#figure-30">Figure 30</a>. The zoom of the `.mzML` - `.featureXML` overlay shows the individual mass traces and the assembly of those in a feature (see <a href="#figure-31">Fig. 31</a>).

(Figure_28)=
|![Opened .mzML in TOPPView](/_images/openms-user-tutorial/metabo/ToppView_1.png)|
|:--:|
|Figure 28: Opened .mzML in TOPPView|

(Figure_29)=
|![Add new layer in TOPPView](/_images/openms-user-tutorial/metabo/ToppView_2.png)|
|:--:|
|Figure 29: Add new layer in TOPPView|

(Figure_30)=
|![Overlay of the .mzML layer with the .featureXML layer](/_images/openms-user-tutorial/metabo/ToppView_3.png)|
|:--:|
|Figure 30: Overlay of the .mzML layer with the .featureXML layer|

(Figure_31)=
|![Zoom of the overlay of the .mzML with the .featureXML layer](/_images/openms-user-tutorial/metabo/ToppView_4.png)|
|:--:|
|Figure 31: Zoom of the overlay of the .mzML with the .featureXML layer. Here the individual isotope traces (blue lines) are assembled into a feature here shown as convex hull (rectangular box).|

To correct precursor m/z values in centroided high-resolution data,, one can add a `HighResPrecursorMassCorrector`node (**Community Nodes**>**OpenMS**>**Mass Correction and Calibration**) before the `FeatureFinderMetabo` node. The following parameter will be used to refine precursor masses by selecting the most intense centroided MS1 peak within a m/z tolerance of 10ppm:

|**parameter**|	**value**|
|:------------|:---------|
|*highest_intensity_peak* → *mz_tolerance* → |100.0|
|*highest_intensity_peak* → *mz_tolerance_unit* → |ppm|

The parameters change the behavior of `HighResPrecursorMassCorrector` as follows:
- **mz_tolerance**: The precursor mass tolerance to find the highest intensity MS1 peak.
- **mz_tolerance_unit**: Unit of precursor mass tolerance.

The workflow can be extended for multi-file analysis, here an `Input Files` node is to be used instead of the `File Importer` node. In front of the `HighResPrecursorMassCorrector`, a `ZipLoopStart` and behind `ZipLoopEnd` has to be used, since `HighResPrecursorMassCorrector` and `FeatureFinderMetabo` will analysis on file to file bases.

To facilitate the collection of features corresponding to the same compound ion across different samples, an alignment of the samples’ feature maps along retention time is often helpful. In addition to local, small-scale elution differences, one can often see constant retention time shifts across large sections between samples. We can use linear transformations to correct for these large scale retention differences. This brings the majority of corresponding compound ions close to each other. Finding the correct corresponding ions is then faster and easier, as we don’t have to search as far around individual features.

(Figure_32)=
|![map alignment example](/_images/openms-user-tutorial/metabo/align.png)|
|:--:|
|Figure 32: The first feature map is used as a reference to which other maps are aligned. The calculated transformation brings corresponding features into close retention time proximity. Linking of these features form a so-called consensus features of a consensus map.|

- After the `ZipLoopEnd` node, add a `MapAlignerPoseClustering` node (**Community Nodes**>**OpenMS**>**Map Alignment**), set its Output Type to featureXML, and adjust the following settings:

|**parameter**|	**value**|
|:------------|:---------|
|*algorithm* → *pairfinder* → *distance_MZ* → *max_difference*|10.0|
|*algorithm* → *pairfinder* → *distance_MZ* → *unit*|ppm|

`MapAlignerPoseClustering` provides an algorithm to align the retention time scales of multiple input files, correcting shifts and distortions between them. Retention time adjustment may be necessary to correct for chromatography differences e.g. before data from multiple LC-MS runs can be combined (feature linking). The alignment algorithm implemented here is the pose clustering algorithm.

The parameters change the behavior of `MapAlignerPoseClustering` as follows:
- **max_difference**: Features that have a larger m/z difference will never be paired.
- **unit**: Unit used for the parameter distance_MZ max_difference, either Da or ppm.

The next step after retention time correction is the grouping of corresponding features in multiple samples. In contrast to the previous alignment, we assume no linear relations of features across samples. The used method is tolerant against local swaps in elution order.

(Figure_33)=
|![feature linking example](/_images/openms-user-tutorial/metabo/link.png)|
|:--:|
|Figure 33: Features A and B correspond to the same analyte. The linking of features between runs (indicated by an arrow) allows comparing feature intensities.|

- After the `MapAlignerPoseClustering` node, add a `FeatureLinkerUnlabeledKD` node (**Community Nodes** > **OpenMS**>**Map Alignment**) and adjust the following settings:

  |**parameter**|**value**|
        |:------------|:--------|
  |*algorithm* → *warp* → *enabled*|false|
  |*algorithm* → *link* → *rt_tol*|30|
  |*algorithm* → *link* → *mz_tol*|10|

  The parameters change the behavior of `FeatureLinkerUnlabeledKD` as follows (similar to the parameters we adjusted for `MapAlignerPoseClustering`):

    - **warp → enabled**: If set to true, feature RTs are warped using LOWESS transformation before linking (reported RTs in results will always be the original RTs).
    - **link → rt_tol**: Width of RT tolerance window (sec).
    - **link → mz_tol**: M/z tolerance.

- After the `FeatureLinkerUnlabeledKD` node, add a **TextExporter** node (**Community Nodes** > **OpenMS** > **File Handling**).
- Add an `Output Folder` node and configure it with an output directory where you want to store the resulting files.
- Run the pipeline and inspect the output.

(Figure_34)=
|![Label-free quantification workflow for metabolites](/_images/openms-user-tutorial/metabo/metabo_part1_with_labels.png)|
|:--:|
|Figure 34: Label-free quantification workflow for metabolites.|

You should find a single, tab-separated file containing the information on where metabolites were found and with which intensities. You can also add `Output Folder` nodes at different stages of the workflow and inspect the intermediate results (e.g., identified metabolite features for each input map). The complete workflow can be seen in <a href="#figure-34">Figure 34</a>. In the following section we will try to identify those metabolites.

<!-- Here the students will have to identify the metabolite of interest and try to find it in the data. -->

The `FeatureLinkerUnlabeledKD` output can be visualized in TOPPView on top of the input and output of the `FeatureFinderMetabo` (see <a href="#figure-35">Fig 35</a>).

(Figure_35)=
|![Label-free quantification workflow for metabolites](/_images/openms-user-tutorial/metabo/ToppView_5.png)|
|:--:|
|Figure 35: Visualization of .consensusXML output over the .mzML and .featureXML ’layer’.|

## Basic metabolite identification
<!-- TODO: Ensure this works with new dataset, at least changing metabolite db to B. subtilis is needed -->
<!-- Might Remove Start -->
At the current state we found several metabolites in the individual maps but so far don’t know what they are. To identify metabolites, OpenMS provides multiple tools, including search by mass: the AccurateMassSearch node searches observed masses against the Human Metabolome Database (HMDB)[^1]<sup>,</sup> [^2]<sup>,</sup> [^3]. We start with the workflow from the previous section (see <a href="#figure-34">Figure 34</a>).

- Add a **FileConverter** node (**Community Nodes** > **OpenMS** > **File Handling**) and connect the output of the FeatureLinkerUnlabeledKD to the incoming port.
- Open the Configure dialog of the **FileConverter** node and select the tab **OutputTypes**. In the drop down list for FileConverter.1.out select **featureXML**.
- Add an **AccurateMassSearch** node (**Community Nodes** > **OpenMS** > **Utilities**) and connect the output of the **FileConverter** node to the first port of the **AccurateMassSearch** node.
- Add four `File Importer` nodes and configure them with the following [files](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Tutorials/Example_Data/Metabolomics/databases/):
    - {path}`Example_Data,Metabolomics,databases,PositiveAdducts.tsv`
      This file specifies the list of adducts that are considered in the positive mode. Each line contains the formula and charge of an adduct separated by a semicolon (e.g. M+H;1+). The mass of the adduct is calculated automatically.
    - {path}`Example_Data,Metabolomics,databases,NegativeAdducts.tsv`
      This file specifies the list of adducts that are considered in the negative mode analogous to the positive mode.
    - {path}`Example_Data,Metabolomics,databases,HMDBMappingFile.tsv`
      This file contains information from a metabolite database in this case from HMDB. It has three (or more) tab-separated columns: mass, formula, and identifier(s). This allows for an efficient search by mass.
    - {path}`Example_Data,Metabolomics,databases,HMDB2StructMapping.tsv`
      This file contains additional information about the identifiers in the mapping file. It has four tab-separated columns that contain the identifier, name, SMILES, and INCHI. These will be included in the result file. The identifiers in this file must match the identifiers in the `HMDBMappingFile.tsv`.
- In the same order as they are given above connect them to the remaining input ports of the **AccurateMassSearch** node.
- Add an `Output Folder` node and connect the first output port of the
  **AccurateMassSearch** node to the `Output Folder` node.

The result of the **AccurateMassSearch** node is in the mzTab format[^4] so you can easily open it in a text editor or import it into Excel or KNIME, which we will do in the next section. The complete workflow from this section is shown in <a href="#figure-36">Figure 36</a>.

(Figure_36)=
|![Label-free quantification and identification workflow for metabolites](/_images/openms-user-tutorial/metabo/metabo_part2.png)|
|:--:|
|Figure 36: Label-free quantification and identification workflow for metabolites.|

### Convert your data into a KNIME table

The result from the TextExporter node as well as the result from the **AccurateMassSearch** node are files while standard KNIME nodes display and process only KNIME tables. To convert these files into KNIME tables we need two different nodes. For the **AccurateMassSearch** results, we use the **MzTabReader** node (**Community Nodes** > **OpenMS** > **Conversion** > **mzTab**) and its **Small Molecule Section** port. For the result of the **TextExporter**, we use the `ConsensusTextReader` (**Community Nodes** > **OpenMS** > **Conversion**).
When executed, both nodes will import the OpenMS files and provide access to the data as KNIME tables. The retention time values are exported as a list using the **MzTabReader** based on the current PSI-Standard. This has to be parsed using the **SplitCollectionColumn**, which outputs a ”Split Value 1” based on the first entry in the rention time list, which has to be renamed to retention time using the **ColumnRename**. You can now combine both tables using the **Joiner** node (**Manipulation** > **Column** > **Split & Combine**) and configure it to match the m/z and retention time values of the respective tables. The full workflow is shown in <a href="#figure-37">Figure 37</a>.

(Figure_37)=
|![Label-free quantification and identification workflow for metabolites that loads the results into KNIME and joins the tables](/_images/openms-user-tutorial/metabo/metabo_part3.png)|
|:--:|
|Figure 37: Label-free quantification and identification workflow for metabolites that loads the results into KNIME and joins the tables.|

### Adduct grouping

Metabolites commonly co-elute as ions with different adducts (e.g., glutathione+H, glutathione+Na) or with charge-neutral modifications (e.g., water loss). Grouping such related ions allows to leverage information across features. For example, a low intensity, single trace feature could still be assigned a charge and adduct due to a matching high-quality feature. Several OpenMS tools, such as **AccurateMassSearch**, can use this information to, for example, narrow down candidates for identification.

For this grouping task, we provide the **MetaboliteAdductDecharger** node. Its method explores the combinatorial space of all adduct combinations in a charge range for optimal explanations. Using defined adduct probabilities, it assigns co-eluting features having suitable mass shifts and charges those adduct combinations which maximize overall ion probabilities.

The tool works natively with featureXML data, allowing the use of reported convex hulls. On such a single-sample level, co-elution settings can be chosen more stringently, as ionization-based adducts should not influence the elution time: Instead, elution differences of related ions should be due to slightly differently estimated times for their feature centroids.

Alternatively, consensusXML data from feature linking can be converted for use, though with less chromatographic information. Here, the elution time averaging for features linked across samples, motivates wider co-elution tolerances.

The two main tool outputs are a consensusXML file with compound groups of related input ions, and a featureXML containing the input file but annotated with inferred adduct information and charges.

Options to respect or replace ion charges or adducts allow for example:
- Heuristic but faster, iterative adduct grouping(**MetaboliteAdductDecharger → MetaboliteFeatureDeconvolution → q_try** set to “feature”) by chaining multiple **MetaboliteAdductDecharger** nodes with growing adduct sets, charge ranges or otherwise relaxed tolerances.
- More specific feature linking (**FeatureLinkerUnlabeledKD → algorithm → ignore_adduct** set to “false”)

(Figure_38)=
|![Metabolite Adduct Decharger adduct grouping workflow](/_images/openms-user-tutorial/metabo/mad.png)|
|:--:|
|Figure 38: Metabolite Adduct Decharger adduct grouping workflow. |

<div class="admonition task">
<p class="admonition-title task-title">**Task**</p>
A modified metabolomics workflow with exemplary MetaboliteAdductDecharger use and parameters is provided in {path}`Workflows,MetaboliteAdductGrouping.knwf`. Run the workflow, inspect tool outputs and compare **AccurateMassSearch** results with and without adduct grouping.
</div>

### Visualizing data

Now that you have your data in KNIME you should try to get a feeling for the capabilities of KNIME.

<div class="admonition task">
<p class="admonition-title task-title">**Task**</p>
Check out the **Molecule Type Cast** node (**Chemistry** > **Translators**) together with subsequent cheminformatics nodes (e.g. **RDKit From Molecule**(**Community Nodes** > **RDKit** > **Converters**)) to render the structural formula contained in the result table.
</div>

<div class="admonition task">
<p class="admonition-title task-title">**Task**</p>
Have a look at the `Column Filter` node to reduce the table to the interesting columns, e.g., only the Ids, chemical formula, and intensities.
</div>

<div class="admonition task">
<p class="admonition-title task-title">**Task**</p>
Try to compute and visualize the m/z and retention time error of the different feature elements (from the input maps) of each consensus feature. Hint: A nicely configured **Math Formula (Multi Column)** node should suffice.
</div>

<!-- Might Remove End -->

## Identifying Metabolites Using Spectral Libraries

Relying solely on accurate mass for metabolite identification can yield ambiguous results. In practice, additional information—such as retention time—is used to refine candidate selection. Beyond MS1-based features, tandem mass spectra (MS2) provide crucial structural insights. In this section, we explore how metabolites can be identified using a library of previously characterized spectra.

### Annotating MS1 Features with Fragment Spectra

Before performing a spectral library search, we must annotate MS1-based features with their corresponding fragment spectra. To achieve this, we introduce the **IDMapper** node between **MapAlignerPoseClusteringKD** and **FeatureLinkerUnlabeledKD**. Since **IDMapper** operates on a per-spectrum basis, we loop over the individual `featureXML` files using a **ZipLoop**. Additionally, as we aim to annotate MS spectra, we also connect the original mzML **FileImporter** node to the loop.

Originally designed for proteomics, **IDMapper** expects protein identifications. To bypass this requirement, we connect an **empty ID file** via a **FileImporter** node.

### Filtering Unannotated MS1 Features

After linking features, we remove MS1-based features that lack fragment spectrum annotations. This is done by connecting a **FileFilter** node to the output of **FeatureLinkerUnlabeledKD** and enabling:

```
id → remove_unannotated_features = true
```

### Generating Consensus Spectra

Next, we consolidate individual MS2 spectra into **consensus spectra** using the **GNPSExport** node, which outputs an `mgf` file containing the consensus spectra. To configure this step:

- Connect the **FileFilter** node output to **GNPSExport**.
- Also, connect the **FileImporter** node containing the original spectra.

### Performing Spectral Library Matching

Now, we are ready to identify metabolites using the **MetaboliteSpectralMatcher** node:

1. Add a **FileImporter** node to load the spectral database.
2. Convert the `mgf` output of **GNPSExport** into `mzML` format using the **FileConverter** node, ensuring:
   ```
   algorithm → merge_spectra = false
   ```
3. Read the resulting `mzTab` file into a KNIME table using the **mzTabReader** node.

This completes the workflow for spectral library-based metabolite identification.

## References

[^1]: D. S. Wishart, D. Tzur, C. Knox, et al., HMDB: the Human Metabolome Database,
Nucleic Acids Res 35(Database issue), D521–6 (Jan 2007), <a href="https://academic.oup.com/nar/article/35/suppl_1/D521/1109186">doi:10.1093/nar/gkl923</a>. 69

[^2]: D. S. Wishart, C. Knox, A. C. Guo, et al., HMDB: a knowledgebase for the human
metabolome, Nucleic Acids Res 37(Database issue), D603–10 (Jan 2009), <a href="https://academic.oup.com/nar/article/37/suppl_1/D603/1011821">doi: 10.1093/nar/gkn810</a>. 69

[^3]: D. S. Wishart, T. Jewison, A. C. Guo, M. Wilson, C. Knox, et al., HMDB 3.0–The
Human Metabolome Database in 2013, Nucleic Acids Res 41(Database issue),D801–7 (Jan 2013), <a href="https://academic.oup.com/nar/article/41/D1/D801/1055560">doi:10.1093/nar/gks1065</a>. 69

[^4]: J. Griss, A. R. Jones, T. Sachsenberg, M. Walzer, L. Gatto, J. Hartler, G. G.
Thallinger, R. M. Salek, C. Steinbeck, N. Neuhauser, J. Cox, S. Neumann, J. Fan,
F. Reisinger, Q.-W. Xu, N. Del Toro, Y. Perez-Riverol, F. Ghali, N. Bandeira, I. Xenarios, O. Kohlbacher, J. A. Vizcaino, and H. Hermjakob, The mzTab Data Exchange Format: communicating MS-based proteomics and metabolomics experimental results to a wider audience, Mol Cell Proteomics (Jun 2014), doi:10.1074/mcp.O113.036681. 69

[^5]: S. Böcker, M. C. Letzel, Z. Lipták, and A. Pervukhin, SIRIUS: Decomposing isotope
patterns for metabolite identification, Bioinformatics 25(2), 218–224 (2009), <a href="https://academic.oup.com/bioinformatics/article/25/2/218/218950">doi:10.1093/bioinformatics/btn603</a>. 75

[^6]: S. Böcker and K. Dührkop, Fragmentation trees reloaded, J. Cheminform. 8(1),
1–26 (2016), <a href="https://jcheminf.biomedcentral.com/articles/10.1186/s13321-016-0116-8">doi:10.1186/s13321-016-0116-8</a>. 75

[^7]: K. Dührkop, H. Shen, M. Meusel, J. Rousu, and S. Böcker, <a href="https://www.pnas.org/doi/abs/10.1073/pnas.1509788112">Searching molecular structure databases with tandem mass spectra using CSI:FingerID</a>, Proc. Natl.
Acad. Sci. 112(41), 12580–12585 (oct 2015), <a href="https://www.pnas.org/doi/full/10.1073/pnas.1509788112">doi:10.1073/pnas.1509788112</a>. 75
