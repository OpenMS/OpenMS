MSStats
=======

## R integration

KNIME provides a large number of nodes for a wide range of statistical analysis, machine learning, data processing, and
visualization. Still, more recent statistical analysis methods, specialized visualizations or cutting edge algorithms
may not be covered in KNIME. In order to expand its capabilities beyond the readily available nodes, external scripting
languages can be integrated. In this tutorial, we primarily use scripts of the powerful statistical computing language R.
Note that this part is considered advanced and might be difficult to follow if you are not familiar with R. In this case
you might skip this part.

**R View (Table)** allows to seamlessly include R scripts into KNIME. We will
demonstrate on a minimal example how such a script is integrated.

<div class="admonition task">
  <p class="admonition-title task-title">**Task**</p>
  <p>
First we need some example data in KNIME, which we will generate using the **Data Generator** node (**IO** > **Other** > **Data Generator**).
You can keep the default settings and execute the node. The table contains four columns, each containing random coordinates and one column
containing a cluster number (Cluster_0 to Cluster_3). Now place a **R View (Table)** node into the workflow and connect
the upper output port of the **Data Generator** node to the input of the **R View (Table)** node. Right-click and
configure the node. If you get an error message like `Execute failed: R_HOME does not contain a folder with name ’bin’.`
or `Execution failed: R Home is invalid.`: please change the R settings in the preferences. To do so open **File** >
**Preferences** > **KNIME** > **R** and enter the path to your R installation (the folder that contains the bin
directory. e.g., {path}`C:,Program Files,R,R-3.4.3`).
  </p>
  <p>
If you get an error message like: ”Execute failed: Could not find Rserve package. Please install it in your R
installation by running ”install.packages(’Rserve’)”.” You may need to run your R binary as administrator (In windows
explorer: right-click ”Run as administrator”) and enter install.packages(’Rserve’) to install the package.
  </p>
  <p>
If R is correctly recognized we can start writing an R script. Consider that we are interested in plotting the first and
second coordinates and color them according to their cluster number. In R this can be done in a single line. In the
**R view (Table)** text editor, enter the following code:
```r
plot(x=knime.in$Universe_0_0, y=knime.in$Universe_0_1, main="Plotting column Universe_0_0 vs. Universe_0_1", col=knime.in$"Cluster Membership")
```
  </p>
  <p>
**Explanation:** The table provided as input to the **R View (Table)** node is available as R **data.frame** with name
`knime.in`. Columns (also listed on the left side of the R View window) can be accessed in the usual R way by first
specifying the `data.frame` name and then the column name (e.g., `knime.in$Universe_0_0`). `plot` is the plotting function
we use to generate the image. We tell it to use the data in column `Universe_0_0` of the dataframe object **knime.in**
(denoted as `knime.in$Universe_0_0`) as x-coordinate and the other column `knime.in$Universe_0_1` as y-coordinate in the
plot. `main` is simply the main title of the plot and `col` the column that is used to determine the color (in this case
it is the `Cluster Membership` column).
  </p>
  <p>
Now press the <kbd>Eval script</kbd> and <kbd>Show plot</kbd> buttons.
  </p>
</div>

```{note}
Note that we needed to put some extra quotes around `Cluster Membership`. If we omit those, R would interpret the column
name only up to the first space `(knime.in$Cluster)` which is not present in the table and leads to an error. Quotes are
regularly needed if column names contain spaces, tabs or other special characters like $ itself.
```

## Using MSstats in a KNIME workflow

The R package `MSstats` can be used for statistical relative quantification of proteins and peptides in mass spectrometry-based proteomics. Supported are label-free as well as labeled experiments in combination with data-dependent, targeted and data independent acquisition. Inputs can be identified and quantified entities (peptides or proteins) and the output is a list of differentially abundant entities, or summaries of their relative abundance. It depends on accurate feature detection, identification
and quantification which can be performed e.g. by an OpenMS workflow. MSstats can be used for data processing & visualization, as well as statistical modeling & inference. Please see [^1] and the [MSstats](http://msstats.org) website for further
information.

### Identification and quantification of the iPRG2015 data with subsequent MSstats analysis

Here, we describe how to use OpenMS and MSstats for the analysis of the ABRF iPRG2015 dataset[^2].

```{note}
Reanalysing the full dataset from scratch would take too long. In the following tutorial, we will focus on just the conversion process and the downstream analysis.
```


#### Dataset

The iPRG (Proteome Informatics Research Group) dataset from the study in 2015, as
described in [^2], aims at evaluating the effect of statistical analysis software on the
accuracy of results on a proteomics label-free quantification experiment. The data is
based on four artificial samples with known composition (background: 200 ng *S. cerevisiae*). These were spiked with different quantities of individual digested proteins,
whose identifiers were masked for the competition as yeast proteins in the provided
database (see <a href="#table-1">Table 1</a>).

<div class="table" id="table-1">

<!-- l. 1198 --><p class="indent">   </p><figure class="float">

<a id="x1-32001r1"></a>
<a id="x1-32002"></a>
<figcaption class="caption"><span class="id">Table&nbsp;1: Samples (background: 200 ng <i>S.&nbsp;cerevisiae</i>) with spiked-in proteins <br> in different
quantities [fmols]</span></figcaption><!-- tex4ht:label?: x1-32001r3  -->
<div class="tabular">
 <table class="tabular" id="TBL-1"><colgroup id="TBL-1-1g"><col id="TBL-1-1"><col id="TBL-1-2"><col id="TBL-1-3"><col id="TBL-1-4"><col id="TBL-1-5"><col id="TBL-1-6"><col id="TBL-1-7"><col id="TBL-1-8"></colgroup><tbody><tr id="TBL-1-1-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-1-1">   </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-1-2">                      </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-1-3">                   </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-1-4">                 </td><td class="td11" style="white-space:nowrap; text-align:right;" colspan="4" id="TBL-1-1-5">         <div class="multicolumn" style="white-space:nowrap; text-align:center;"><span class="rm-lmr-10x-x-109">Samples</span></div>
</td></tr><tr id="TBL-1-2-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-2-1">    </td><td class="td11" style="white-space:nowrap; text-align:center;" id="TBL-1-2-2">        <div class="multicolumn" style="white-space:nowrap; text-align:center;"><span class="rm-lmr-10x-x-109">Name</span></div>         </td><td class="td11" style="white-space:nowrap; text-align:center;" id="TBL-1-2-3">       <div class="multicolumn" style="white-space:nowrap; text-align:center;"><span class="rm-lmr-10x-x-109">Origin</span></div>       </td><td class="td11" style="white-space:nowrap; text-align:center;" id="TBL-1-2-4"> <div class="multicolumn" style="white-space:nowrap; text-align:center;"><span class="rm-lmr-10x-x-109">Molecular Weight</span></div> </td><td class="td11" style="white-space:nowrap; text-align:center;" id="TBL-1-2-5">  <div class="multicolumn" style="white-space:nowrap; text-align:center;"><span class="rm-lmr-10x-x-109">1</span></div>  </td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-2-6">   <span class="rm-lmr-10x-x-109">2  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-2-7">  <span class="rm-lmr-10x-x-109">3  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-2-8">   <span class="rm-lmr-10x-x-109">4  </span></td>
</tr><tr id="TBL-1-3-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-3-1"> <span class="rm-lmr-10x-x-109">A  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-3-2"> <span class="rm-lmr-10x-x-109">Ovalbumin                  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-3-3"> <i>Egg White</i>             </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-3-4"> <span class="rm-lmr-10x-x-109">45 KD                </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-3-5"> <span class="rm-lmr-10x-x-109">65  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-3-6">  <span class="rm-lmr-10x-x-109">55  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-3-7"> <span class="rm-lmr-10x-x-109">15  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-3-8">   <span class="rm-lmr-10x-x-109">2  </span></td>
</tr><tr id="TBL-1-4-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-4-1"> <span class="rm-lmr-10x-x-109">B  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-4-2"> <span class="rm-lmr-10x-x-109">Myoglobin                   </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-4-3"> <i>Equine Heart</i>          </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-4-4"> <span class="rm-lmr-10x-x-109">17 KD                </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-4-5"> <span class="rm-lmr-10x-x-109">55  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-4-6">  <span class="rm-lmr-10x-x-109">15  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-4-7">  <span class="rm-lmr-10x-x-109">2  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-4-8">  <span class="rm-lmr-10x-x-109">65  </span></td>
</tr><tr id="TBL-1-5-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-5-1"> <span class="rm-lmr-10x-x-109">C  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-5-2"> <span class="rm-lmr-10x-x-109">Phosphorylase b           </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-5-3"> <i>Rabbit Muscle</i>         </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-5-4"> <span class="rm-lmr-10x-x-109">97 KD                </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-5-5"> <span class="rm-lmr-10x-x-109">15  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-5-6">   <span class="rm-lmr-10x-x-109">2  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-5-7"> <span class="rm-lmr-10x-x-109">65  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-5-8">  <span class="rm-lmr-10x-x-109">55  </span></td>
</tr><tr id="TBL-1-6-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-6-1"> <span class="rm-lmr-10x-x-109">D  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-6-2"> <span class="rm-lmr-10x-x-109">Beta-Glactosidase         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-6-3"> <i>Escherichia Coli</i>      </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-6-4"> <span class="rm-lmr-10x-x-109">116 KD               </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-6-5">  <span class="rm-lmr-10x-x-109">2  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-6-6">  <span class="rm-lmr-10x-x-109">65  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-6-7"> <span class="rm-lmr-10x-x-109">55  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-6-8">  <span class="rm-lmr-10x-x-109">15  </span></td>
</tr><tr id="TBL-1-7-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-7-1"> <span class="rm-lmr-10x-x-109">E  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-7-2"> <span class="rm-lmr-10x-x-109">Bovine Serum Albumin  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-7-3"> <i>Bovine Serum</i>         </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-7-4"> <span class="rm-lmr-10x-x-109">66 KD                </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-7-5"> <span class="rm-lmr-10x-x-109">11  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-7-6"> <span class="rm-lmr-10x-x-109">0.6  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-7-7"> <span class="rm-lmr-10x-x-109">10  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-7-8"> <span class="rm-lmr-10x-x-109">500  </span></td>
</tr><tr id="TBL-1-8-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-8-1"> <span class="rm-lmr-10x-x-109">F  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-8-2"> <span class="rm-lmr-10x-x-109">Carbonic Anhydrase     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-8-3"> <i>Bovine Erythrocytes</i>  </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-1-8-4"> <span class="rm-lmr-10x-x-109">29 KD                </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-8-5"> <span class="rm-lmr-10x-x-109">10  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-8-6"> <span class="rm-lmr-10x-x-109">500  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-8-7"> <span class="rm-lmr-10x-x-109">11  </span></td><td class="td11" style="white-space:nowrap; text-align:right;" id="TBL-1-8-8"> <span class="rm-lmr-10x-x-109">0.6  </span></td>

</tr></tbody></table></div>

</figure>
</div>

#### Identification and quantification

|![KNIME data analysis of iPRG LFQ data](/_images/openms-user-tutorial/labelfree/iPRG_lfq.png)|
|:--:|
|Figure 18: KNIME data analysis of iPRG LFQ data.|

The iPRG LFQ workflow (<a href="#figure-18">Fig. 18</a>) consists of an identification and a quantification part. The identification is achieved by searching the computationally calculated MS2 spectra from a sequence database (`File Importer` node, here with the given database from iPRG:
{path}`ExampleData,iPRG2015,database,iPRG2015targetdecoynocontaminants.fasta`
against the MS2 from the original data (`Input Files` node with all mzMLs following {path}`ExampleData,iPRG2015,datasets,JD06232014sample*.mzML` using the `CometAdapter`.

```{note}
If you want to reproduce the results at home, you have to download the iPRG data in mzML format and perform peak picking on it or convert and pick the raw data with `msconvert`.
```

Afterwards, the results are scored using the **FalseDiscoveryRate** node and filtered to obtain only unique peptides (`IDFilter`) since `MSstats` does not support shared peptides, yet. The quantification is achieved by using the `FeatureFinderCentroided` node, which performs the feature detection on the samples (maps). In the end the quantification results are combined with the filtered identification results (`IDMapper`). In addition, a linear retention time alignment is performed (`MapAlignerPoseClustering`), followed by the feature linking process (**FeatureLinkerUnlabledQT**). The **ConsensusMapNormalizer**s is used to normalize the intensities via robust regression over a set of maps and the `IDConflictResolver` assures that only one identification (best score) is associated with a feature. The output of this workflow is a consensusXML file, which can now be converted using the `MSStatsConverter` (see <a href="#conversion-and-downstream-analysis">Conversion and downstream analysis</a> section).

#### Experimental design

As mentioned before, the downstream analysis can be performed using `MSstats`. In this case, an experimental design has to be specified for the OpenMS workflow. The structure of the experimental design used in OpenMS in case of the iPRG dataset is specified in <a href="#table-2">Table 2</a>.

<div class="table" id="table-2">


<!-- l. 1238 --><p class="indent">   </p><figure class="float">


<a id="x1-34001r2"></a>
<a id="x1-34002"></a>
<figcaption class="caption"><span class="id">Table&nbsp;2: OpenMS Experimental design for the iPRG2015 dataset.</span></figcaption><!-- tex4ht:label?: x1-34001r3  -->
<div class="tabular"> <table class="tabular" id="TBL-2"><colgroup id="TBL-2-1g"><col id="TBL-2-1"><col id="TBL-2-2"><col id="TBL-2-3"><col id="TBL-2-4"><col id="TBL-2-5"></colgroup><tbody><tr id="TBL-2-1-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-1-1"> <span class="rm-lmr-10x-x-109">Fraction_Group  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-1-2"> <span class="rm-lmr-10x-x-109">Fraction                </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-1-3"> <span class="rm-lmr-10x-x-109">Spectra_Filepath        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-1-4"> <span class="rm-lmr-10x-x-109">Label  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-1-5"> <span class="rm-lmr-10x-x-109">Sample  </span></td>
</tr><tr id="TBL-2-2-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-2-1"> <span class="rm-lmr-10x-x-109">1                      </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-2-2"> <span class="rm-lmr-10x-x-109">1                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-2-3"> <span class="rm-lmr-10x-x-109">Sample1-A                </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-2-4"> <span class="rm-lmr-10x-x-109">1       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-2-5"> <span class="rm-lmr-10x-x-109">1          </span></td>
</tr><tr id="TBL-2-3-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-3-1"> <span class="rm-lmr-10x-x-109">2                      </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-3-2"> <span class="rm-lmr-10x-x-109">1                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-3-3"> <span class="rm-lmr-10x-x-109">Sample1-B                </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-3-4"> <span class="rm-lmr-10x-x-109">1       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-3-5"> <span class="rm-lmr-10x-x-109">2          </span></td>
</tr><tr id="TBL-2-4-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-4-1"> <span class="rm-lmr-10x-x-109">3                      </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-4-2"> <span class="rm-lmr-10x-x-109">1                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-4-3"> <span class="rm-lmr-10x-x-109">Sample1-C                </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-4-4"> <span class="rm-lmr-10x-x-109">1       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-4-5"> <span class="rm-lmr-10x-x-109">3          </span></td>
</tr><tr id="TBL-2-5-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-5-1"> <span class="rm-lmr-10x-x-109">4                      </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-5-2"> <span class="rm-lmr-10x-x-109">1                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-5-3"> <span class="rm-lmr-10x-x-109">Sample2-A                </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-5-4"> <span class="rm-lmr-10x-x-109">1       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-5-5"> <span class="rm-lmr-10x-x-109">4          </span></td>
</tr><tr id="TBL-2-6-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-6-1"> <span class="rm-lmr-10x-x-109">5                      </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-6-2"> <span class="rm-lmr-10x-x-109">1                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-6-3"> <span class="rm-lmr-10x-x-109">Sample2-B                </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-6-4"> <span class="rm-lmr-10x-x-109">1       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-6-5"> <span class="rm-lmr-10x-x-109">5          </span></td>
</tr><tr id="TBL-2-7-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-7-1"> <span class="rm-lmr-10x-x-109">6                      </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-7-2"> <span class="rm-lmr-10x-x-109">1                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-7-3"> <span class="rm-lmr-10x-x-109">Sample2-C                </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-7-4"> <span class="rm-lmr-10x-x-109">1       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-7-5"> <span class="rm-lmr-10x-x-109">6          </span></td>
</tr><tr id="TBL-2-8-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-8-1"> <span class="rm-lmr-10x-x-109">7                      </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-8-2"> <span class="rm-lmr-10x-x-109">1                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-8-3"> <span class="rm-lmr-10x-x-109">Sample3-A                </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-8-4"> <span class="rm-lmr-10x-x-109">1       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-8-5"> <span class="rm-lmr-10x-x-109">7          </span></td>
</tr><tr id="TBL-2-9-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-9-1"> <span class="rm-lmr-10x-x-109">8                      </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-9-2"> <span class="rm-lmr-10x-x-109">1                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-9-3"> <span class="rm-lmr-10x-x-109">Sample3-B                </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-9-4"> <span class="rm-lmr-10x-x-109">1       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-9-5"> <span class="rm-lmr-10x-x-109">8          </span></td>
</tr><tr id="TBL-2-10-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-10-1"> <span class="rm-lmr-10x-x-109">9                      </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-10-2"> <span class="rm-lmr-10x-x-109">1                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-10-3"> <span class="rm-lmr-10x-x-109">Sample3-C                </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-10-4"> <span class="rm-lmr-10x-x-109">1       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-10-5"> <span class="rm-lmr-10x-x-109">9          </span></td>
</tr><tr id="TBL-2-11-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-11-1"> <span class="rm-lmr-10x-x-109">10                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-11-2"> <span class="rm-lmr-10x-x-109">1                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-11-3"> <span class="rm-lmr-10x-x-109">Sample4-A                </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-11-4"> <span class="rm-lmr-10x-x-109">1       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-11-5"> <span class="rm-lmr-10x-x-109">10        </span></td>
</tr><tr id="TBL-2-12-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-12-1"> <span class="rm-lmr-10x-x-109">11                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-12-2"> <span class="rm-lmr-10x-x-109">1                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-12-3"> <span class="rm-lmr-10x-x-109">Sample4-B                </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-12-4"> <span class="rm-lmr-10x-x-109">1       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-12-5"> <span class="rm-lmr-10x-x-109">11        </span></td>
</tr><tr id="TBL-2-13-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-13-1"> <span class="rm-lmr-10x-x-109">12                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-13-2"> <span class="rm-lmr-10x-x-109">1                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-13-3"> <span class="rm-lmr-10x-x-109">Sample4-C                </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-13-4"> <span class="rm-lmr-10x-x-109">1       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-13-5"> <span class="rm-lmr-10x-x-109">12        </span></td>
</tr><tr id="TBL-2-14-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-14-1">                </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-14-2">                  </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-14-3">                    </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-14-4">       </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-14-5">        </td>
</tr><tr id="TBL-2-15-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-15-1">                </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-15-2">                  </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-15-3">                    </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-15-4">       </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-15-5">        </td>
</tr><tr id="TBL-2-16-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-16-1"> <span class="rm-lmr-10x-x-109">Sample              </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-16-2"> <span class="rm-lmr-10x-x-109">MSstats_Condition  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-16-3"> <span class="rm-lmr-10x-x-109">MSstats_BioReplicate  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-16-4">       </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-16-5">        </td>
</tr><tr id="TBL-2-17-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-17-1"> <span class="rm-lmr-10x-x-109">1                      </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-17-2"> <span class="rm-lmr-10x-x-109">1                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-17-3"> <span class="rm-lmr-10x-x-109">1                             </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-17-4">       </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-17-5">        </td>
</tr><tr id="TBL-2-18-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-18-1"> <span class="rm-lmr-10x-x-109">2                      </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-18-2"> <span class="rm-lmr-10x-x-109">1                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-18-3"> <span class="rm-lmr-10x-x-109">2                             </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-18-4">       </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-18-5">        </td>
</tr><tr id="TBL-2-19-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-19-1"> <span class="rm-lmr-10x-x-109">3                      </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-19-2"> <span class="rm-lmr-10x-x-109">1                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-19-3"> <span class="rm-lmr-10x-x-109">3                             </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-19-4">       </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-19-5">        </td>
</tr><tr id="TBL-2-20-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-20-1"> <span class="rm-lmr-10x-x-109">4                      </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-20-2"> <span class="rm-lmr-10x-x-109">2                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-20-3"> <span class="rm-lmr-10x-x-109">4                             </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-20-4">       </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-20-5">        </td>
</tr><tr id="TBL-2-21-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-21-1"> <span class="rm-lmr-10x-x-109">5                      </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-21-2"> <span class="rm-lmr-10x-x-109">2                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-21-3"> <span class="rm-lmr-10x-x-109">5                             </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-21-4">       </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-21-5">        </td>
</tr><tr id="TBL-2-22-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-22-1"> <span class="rm-lmr-10x-x-109">6                      </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-22-2"> <span class="rm-lmr-10x-x-109">2                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-22-3"> <span class="rm-lmr-10x-x-109">6                             </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-22-4">       </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-22-5">        </td>
</tr><tr id="TBL-2-23-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-23-1"> <span class="rm-lmr-10x-x-109">7                      </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-23-2"> <span class="rm-lmr-10x-x-109">3                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-23-3"> <span class="rm-lmr-10x-x-109">7                             </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-23-4">       </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-23-5">        </td>
</tr><tr id="TBL-2-24-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-24-1"> <span class="rm-lmr-10x-x-109">8                      </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-24-2"> <span class="rm-lmr-10x-x-109">3                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-24-3"> <span class="rm-lmr-10x-x-109">8                             </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-24-4">       </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-24-5">        </td>
</tr><tr id="TBL-2-25-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-25-1"> <span class="rm-lmr-10x-x-109">9                      </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-25-2"> <span class="rm-lmr-10x-x-109">3                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-25-3"> <span class="rm-lmr-10x-x-109">9                             </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-25-4">       </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-25-5">        </td>
</tr><tr id="TBL-2-26-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-26-1"> <span class="rm-lmr-10x-x-109">10                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-26-2"> <span class="rm-lmr-10x-x-109">4                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-26-3"> <span class="rm-lmr-10x-x-109">10                            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-26-4">       </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-26-5">        </td>
</tr><tr id="TBL-2-27-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-27-1"> <span class="rm-lmr-10x-x-109">11                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-27-2"> <span class="rm-lmr-10x-x-109">4                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-27-3"> <span class="rm-lmr-10x-x-109">11                            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-27-4">       </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-27-5">        </td>
</tr><tr id="TBL-2-28-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-28-1"> <span class="rm-lmr-10x-x-109">12                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-28-2"> <span class="rm-lmr-10x-x-109">4                         </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-28-3"> <span class="rm-lmr-10x-x-109">12                            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-28-4">       </td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-2-28-5">        </td></tr></tbody></table></div>


   </figure>
   </div>

An explanation of the variables can be found in <a href="#table-3">Table 3</a>.

<div class="table">


<!-- l. 1275 --><p class="indent">   </p><figure class="float">


<a id="x1-34003r3"></a>
<a id="x1-34004"></a>
<figcaption class="caption"><span class="id">Table&nbsp;3: Explanation of the column of the experimental design table</span></figcaption><!-- tex4ht:label?: x1-34003r3  -->

 <table class="tabular" id="TBL-5"><colgroup id="TBL-5-1g"><col id="TBL-5-1"><col id="TBL-5-2"></colgroup><tbody><tr id="TBL-5-1-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-5-1-1"> <strong>variables</strong>        </td><td class="td11" style="white-space:normal; text-align:left;" id="TBL-5-1-2"> <!-- l. 1293 --><p class="noindent"><strong>value</strong>                                                                         </p></td>
</tr><tr id="TBL-5-2-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-5-2-1"> <i>Fraction_Group</i>  </td><td class="td11" style="white-space:normal; text-align:left;" id="TBL-5-2-2"> <!-- l. 1293 --><p class="noindent"><span class="rm-lmr-10x-x-109">Index used to group fractions and source files.</span>                                </p></td>
</tr><tr id="TBL-5-3-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-5-3-1"> <i>Fraction</i>           </td><td class="td11" style="white-space:normal; text-align:left;" id="TBL-5-3-2"> <!-- l. 1293 --><p class="noindent"><span class="rm-lmr-10x-x-109">1st, 2nd, .., fraction. Note: All runs must have the same number of
  </span><span class="rm-lmr-10x-x-109">fractions.</span>                                                                                </p></td>
</tr><tr id="TBL-5-4-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-5-4-1"> <i>Spectra_Filepath</i>  </td><td class="td11" style="white-space:normal; text-align:left;" id="TBL-5-4-2"> <!-- l. 1293 --><p class="noindent"><span class="rm-lmr-10x-x-109">Path to mzML files</span>                                                                   </p></td>
</tr><tr id="TBL-5-5-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-5-5-1"> <i>Label</i>               </td><td class="td11" style="white-space:normal; text-align:left;" id="TBL-5-5-2"> <!-- l. 1293 --><p class="noindent"><span class="rm-lmr-10x-x-109">label-free: always 1</span>                                                                   </p></td>
</tr><tr id="TBL-5-6-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-5-6-1"> <i></i>               </td><td class="td11" style="white-space:normal; text-align:left;" id="TBL-5-6-2"> <!-- l. 1293 --><p class="noindent"><span class="rm-lmr-10x-x-109">TMT6Plex: 1...6</span>                                                                      </p></td>
</tr><tr id="TBL-5-7-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-5-7-1"> <i></i>               </td><td class="td11" style="white-space:normal; text-align:left;" id="TBL-5-7-2"> <!-- l. 1293 --><p class="noindent"><span class="rm-lmr-10x-x-109">SILAC with light and heavy: 1..2</span>                                                 </p></td>
</tr><tr id="TBL-5-8-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-5-8-1"> <i>Sample</i>             </td><td class="td11" style="white-space:normal; text-align:left;" id="TBL-5-8-2"> <!-- l. 1293 --><p class="noindent"><span class="rm-lmr-10x-x-109">Index of sample measured in the specified label X, in fraction Y of
  </span><span class="rm-lmr-10x-x-109">fraction group Z.</span>                                                                      </p></td>
</tr><tr id="TBL-5-9-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-5-9-1"> <i>Conditions</i>        </td><td class="td11" style="white-space:normal; text-align:left;" id="TBL-5-9-2"> <!-- l. 1293 --><p class="noindent"><span class="rm-lmr-10x-x-109">Further specification of different conditions (e.g. MSstats_Condition;
  </span><span class="rm-lmr-10x-x-109">MSstats_BioReplicate)</span>                                                              </p></td>

</tr></tbody></table>


   </figure>
   </div>

The conditions are highly dependent on the type of experiment and on which kind of analysis you want to perform. For the `MSstats` analysis the information which sample belongs to which condition and if there are biological replicates are mandatory. This can be specified in further condition columns as explained in Table 3. For a detailed description of the MSstats-specific terminology, see their documentation e.g. in the R vignette.

#### Conversion and downstream analysis

Conversion of the OpenMS-internal consensusXML format (which is an aggregation of quantified and possibly identified features across several MS-maps) to a table (in MSstats-conformant CSV format) is very easy. First, create a new KNIME workflow. Then, run the `MSStatsConverter` node with a consensusXML and the manually created (e.g. in Excel) experimental design as inputs (loaded via `File Importer` nodes). The first input can be found in:

{path}`ExampleData,iPRG2015,openmsLFQResults,iPRGlfq.consensusXML`

This file was generated by using the {path}`Workflows,openmsLFQiPRG2015.knwf` workflow (seen in <a href="#figure-18">Fig. 18</a>). The second input is specified in:

{path}`ExampleData,iPRG2015,experimentaldesign.tsv`

Adjust the parameters in the config dialog of the converter to match the given experimental design file and to use a simple summing for peptides that elute in multiple features (with the same charge state, i.e. m/z value).

|**parameter**|**value**|
|------------:|:--------|
|*msstats_bioreplicate*|MSstats_Bioreplicate|
|*msstats_condition*|MSstats_Condition|
|*labeled_reference_peptides*|false|
|*retention_time_summarization_method (advanced)*|sum|

The downstream analysis of the peptide ions with `MSstats` is performed in several steps. These steps are reflected by several KNIME R nodes, which consume the output of `MSStatsConverter`. The outline of the workflow is shown in <a href="#figure-19">Figure 19</a>.

|![MSstats analysis using KNIME](/_images/openms-user-tutorial/labelfree/MSstats.png)|
|:--:|
|Figure 19: MSstats analysis using KNIME. The individual steps (Preprocessing, Group Comparisons, Result Data Renaming, and Export) are split among several consecutive nodes.|

We load the file resulting from `MSStatsConverter` either by saving it with an `Output File` node and reloading it with the `File Reader`.  Alternatively, for advanced users, you can use a URI Port to Variable node and use the variable in the File Reader config dialog (**V** button - located on the right of the **Browse** button) to read from the temporary file.

##### Preprocessing

The first node (`Table to R`) loads `MSstats` as well as the data from the previous KNIME node and performs a preprocessing step on the input data. The following inline R script ( needs to be pasted into the config dialog of the node):

```r
library(MSstats)
data <- knime.in
quant <- OpenMStoMSstatsFormat(data, removeProtein_with1Feature = FALSE)
```

The inline R script allows further preparation of the data produced by `MSStatsConverter` before the actual analysis is performed. In this example, the lines with proteins, which were identified with only one feature, were retained. Alternatively they could be removed.
In the same node, most importantly, the following line transforms the data into a format that is understood by `MSstats`:

```r
processed.quant <- dataProcess(quant, censoredInt = 'NA')
```
Here, `dataProcess` is one of the most important functions that the R package provides. The function performs the following steps:

1. Logarithm transformation of the intensities
2. Normalization
3. Feature selection
4. Missing value imputation
5. Run-level summarization

In this example, we just state that missing intensity values are represented by the `NA` string.

##### Group Comparison

The goal of the analysis is the determination of differentially-expressed proteins among the different conditions C1-C4. We can specify the comparisons that we want to make in a *comparison* matrix. For this, let’s consider the following example:

![comparison matrix](/_images/openms-user-tutorial/labelfree/handout-clean129x.svg)

This matrix has the following properties:

- The number of rows equals the number of comparisons that we want to perform, the number of columns equals the number of conditions (here, column 1 refers to C1, column 2 to C2 and so forth).
- The entries of each row consist of exactly one 1 and one -1, the others must be 0.
- The condition with the entry `1` constitutes the enumerator of the log2 fold-change. The one with entry `-1` denotes the denominator. Hence, the first row states that we want calculate:
```{math}
 \begin{equation} \log_2 \frac{C_{2}}{C_{1}} \end{equation}
```
We can generate such a matrix in R using the following code snippet in (for example) a new **R to R** node that takes over the R workspace from the previous node with all its variables:

```r
comparison1<-matrix(c(-1,1,0,0),nrow=1)   
comparison2<-matrix(c(-1,0,1,0),nrow=1)

comparison3<-matrix(c(-1,0,0,1),nrow=1)  
comparison4<-matrix(c(0,-1,1,0),nrow=1)

comparison5<-matrix(c(0,-1,0,1),nrow=1)  
comparison6<-matrix(c(0,0,-1,1),nrow=1)

comparison <- rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6)
row.names(comparison)<-c("C2-C1","C3-C1","C4-C1","C3-C2","C4-C2","C4-C3")
```

Here, we assemble each row in turn, concatenate them at the end, and provide row names for labeling the rows with the respective condition.
In MSstats, the group comparison is then performed with the following line:

```r
test.MSstats <- groupComparison(contrast.matrix=comparison, data=processed.quant)
```
No more parameters need to be set for performing the comparison.

##### Result processing

In a next R to R node, the results are being processed. The following code snippet will rename the spiked-in proteins to A,B,C,D,E, and F and remove the names of other proteins, which will be beneficial for the subsequent visualization, as for example performed in <a href="#figure-20">Figure 20</a>:

```r
  test.MSstats.cr <- test.MSstats$ComparisonResult   


  # Rename spiked ins to A,B,C....  
  pnames <- c("A", "B", "C", "D", "E", "F")

  names(pnames) <- c(  
  "sp|P44015|VAC2_YEAST",  
  "sp|P55752|ISCB_YEAST",

  "sp|P44374|SFG2_YEAST",  
  "sp|P44983|UTR6_YEAST",  
  "sp|P44683|PGA4_YEAST",

  "sp|P55249|ZRT4_YEAST"  
  )  

  test.MSstats.cr.spikedins <- bind_rows(

  test.MSstats.cr[grep("P44015", test.MSstats.cr$Protein),],

  test.MSstats.cr[grep("P55752", test.MSstats.cr$Protein),],

  test.MSstats.cr[grep("P44374", test.MSstats.cr$Protein),],

  test.MSstats.cr[grep("P44683", test.MSstats.cr$Protein),],

  test.MSstats.cr[grep("P44983", test.MSstats.cr$Protein),],

  test.MSstats.cr[grep("P55249", test.MSstats.cr$Protein),]  
  )  
  # Rename Proteins

  test.MSstats.cr.spikedins$Protein <- sapply(test.MSstats.cr.spikedins$Protein, function(x) {pnames[as.character(x)]})



  test.MSstats.cr$Protein <- sapply(test.MSstats.cr$Protein, function(x) {


    x <- as.character(x)  

    if (x %in% names(pnames)) {


      return(pnames[as.character(x)])  
      } else {  
      return("")

    }
  })

```

##### Export

The last four nodes, each connected and making use of the same workspace from the last node, will export the results to a textual representation and volcano plots for further inspection. Firstly, quality control can be performed with the following snippet:

```r
qcplot <- dataProcessPlots(processed.quant, type="QCplot",   
        ylimDown=0,

which.Protein = 'allonly',
width=7, height=7, address=F)
```

The code for this snippet is embedded in the first output node of the workflow. The resulting boxplots show the log2 intensity distribution across the MS runs.
The second node is an **R View (Workspace)** node that returns a Volcano plot which displays differentially expressed proteins between conditions C2 vs. C1. The plot is described in more detail in the following Result section. This is how you generate it:

```r
groupComparisonPlots(data=test.MSstats.cr, type="VolcanoPlot",

                    width=12, height=12,dot.size = 2,ylimUp = 7,

                    which.Comparison = "C2-C1",
                    address=F)
```
The last two nodes export the `MSstats` results as a KNIME table for potential further analysis or for writing it to a (e.g. csv) file. Note that you could also write output inside the Rscript if you are familiar with it. Use the following for an **R to Table** node exporting all results:

```r
knime.out <- test.MSstats.cr
```
And this for an **R to Table** node exporting only results for the spike-ins:

```r
knime.out <- test.MSstats.cr.spikedins
```

#### Result

An excerpt of the main result of the group comparison can be seen in <a href="#figure-20">Figure 20</a>.

|![Volcano plots c2_c1](/_images/openms-user-tutorial/labelfree/c2_c1-.png) ![Volcano plots c3_c2](/_images/openms-user-tutorial/labelfree/c3_c2-.png)|
|:--:|
|Figure 20: Volcano plots produced by the Group Comparison in MSstats The dotted line indicates an adjusted p-value threshold|

The Volcano plots show differently expressed spiked-in proteins. In the left plot, which shows the fold-change C2-C1, we can see the proteins D and F (`sp|P44983|UTR6_YEAST` and `sp|P55249|ZRT4_YEAST`) are significantly over-expressed in C2, while the proteins B,C, and E (`sp|P55752|ISCB_YEAST`, `sp|P55752|ISCB_YEAST`, and `sp|P44683|PGA4_YEAST`) are under-expressed. In the right plot, which shows the fold-change ratio of C3 vs. C2, we can see the proteins E and C (`sp|P44683|PGA4_YEAST` and `sp|P44374|SFG2_YEAST`) over-expressed and the proteins A and F (`sp|P44015|VAC2_YEAST` and `sp|P55249|ZRT4_YEAST`) under-expressed. The plots also show further differentially-expressed proteins, which do not belong to the spiked-in proteins.

The full analysis workflow can be found under:
{path}`Workflows,MSstatsstatPostProcessingiPRG2015.knwf`


### Isobaric analysis workflow

In the last chapters, we identified and quantified peptides in a label-free experiment.

In this section, we would like to introduce a possible workflow for the analysis of isobaric data.
Let’s have a look at the workflow (see <a href="#figure-23">Fig 23</a>).

(Figure_23)=
|![Workflow for the analysis of isobaric data](/_images/openms-user-tutorial/isobaric/isobaric_inference_wf.png)|
|:--:|
|Figure 23: Workflow for the analysis of isobaric data|

The full analysis workflow can be found here:
{path}`Workflows,IdentificationquantificationisobaricinferenceepifanyMSstatsTMT`

The workflow has four input nodes. The first for the experimental design to allow for MSstatsTMT compatible export (`MSStatsConverter`). The second for the `.mzML` files with the centroided spectra from the isobaric labeling experiment and the third one for the `.fasta` database used for identification. The last one allows to specify an output path for the plots generated by the R View, which runs MSstatsTMT (I). The quantification (A) is performed using the **IsobaricAnalzyer**. The tool is able to extract and normalize quantitative information from TMT and iTRAQ data. The values can be assessed from centroided MS2 or MS3 spectra (if available). Isotope correction is performed based on the specified correction matrix (as provided by the manufacturer). The identification (C) is applied as known from the previous chapters by using database search and a target-decoy database.

To reduce the complexity of the data for later inference the q-value estimation and FDR filtering is performed on PSM level for each file individually (B). Afterwards the identification (PSM) and quantiative information is combined using the `IDMapper`. After the processing of all available files, the intermediate results are aggregated (**FileMerger** - D). All PSM results are used for score estimation and protein inference (**Epifany**) (E). For detailed information about protein inference please see Chaper 4. Then, decoys are removed and the inference results are filtered via a protein group FDR. Peptide level results can be exported via **MzTabExporter** (F), protein level results can be obtained via the **ProteinQuantifier** (G) or the results can exported (`MSStatsConverter` - H) and further processed with the following R pipeline to allow for downstream processing using `MSstatsTMT`.

Please import the workflow from {path}`Workflows,IdentificationquantificationisobaricinferenceepifanyMSstatsTMT` into KNIME via the menu entry **File** > **Import KNIME workflow** > **Select file** and double-click the imported workflow in order to open it. Before you can execute the workflow, you have to correct the locations of the files in the `Input Files` nodes (don’t forget the one for the FASTA database inside the “ID” meta node). Try and run your workflow by executing all nodes at once.

#### Excursion MSstatsTMT

The R package `MSstatsTMT` can be used for protein significance analysis in shotgun mass spectrometry-based proteomic experiments with tandem mass tag (TMT) labeling. `MSstatsTMT` provides functionality for two types of analysis & their visualization: Protein summarization based on peptide quantification and Model-based group comparison to detect significant changes in abundance. It depends on accurate feature detection, identification and quantification which can be performed e.g. by an OpenMS workflow.

In general, `MSstatsTMT` can be used for data processing & visualization, as well as statistical modeling. Please see [^3] and the [MSstats](http://msstats.org/msstatstmt/) website for further information.

There is also an [online lecture](https://youtu.be/3CDnrQxGLbA) and tutorial for `MSstatsTMT` from the May Institute Workshop 2020.

#### Dataset and experimental design

We are using the MSV000084264 ground truth dataset, which consists of TMT10plex controlled mixes of different concentrated UPS1 peptides spiked into SILAC HeLa peptides measured in a dilution series https://www.omicsdi.org/dataset/massive/MSV000084264. <a href="#figure-24">Figure 24</a> shows the experimental design. In this experiment, 5 different TMT10plex mixtures – different labeling strategies – were analysed. These were measured in triplicates represented by the 15 MS runs (3 runs each). The example data, database and experimental design to run the workflow can be found [here](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Tutorials/Data/isobaric_MSV000084264/).

(Figure_24)=
|![Experimental Design](/_images/openms-user-tutorial/isobaric/isobaric_experimental_design.jpeg)|
|:--:|
|Figure 24: Experimental Design|

The experimental design in table format allows for `MSstatsTMT` compatible export. The design is represented by two tables. The first one <a href="#table-4">4</a> represents the overall structure of the experiment in terms of samples, fractions, labels and fraction groups. The second one <a href="#table-5">5</a> adds to the first by specifying specific conditions, biological replicates as well as mixtures and label for each channel. For additional information about the experimental design please see <a href="#table-3">Table 3</a>.

<div class="table" id="table-4">


<!-- l. 1737 --><p class="indent">   </p><figure class="float">


<a id="x1-50003r4"></a>
<a id="x1-50004"></a>
<figcaption class="caption"><span class="id">Table&nbsp;4: </span><span class="content">Experimental Design 1
</span></figcaption><!-- tex4ht:label?: x1-50003r5  -->
 <table class="tabular" id="TBL-7"><colgroup id="TBL-7-1g"><col id="TBL-7-1"><col id="TBL-7-2"><col id="TBL-7-3"><col id="TBL-7-4"><col id="TBL-7-5"></colgroup><tbody><tr id="TBL-7-1-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-1-1">  <span class="rm-lmr-6">Spectra_Filepath                                                                   </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-1-2">  <span class="rm-lmr-6">Fraction  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-1-3">  <span class="rm-lmr-6">Label  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-1-4">  <span class="rm-lmr-6">Fraction_Group  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-1-5">  <span class="rm-lmr-6">Sample  </span></td>
</tr><tr id="TBL-7-2-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-2-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_01.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-2-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-2-3">  <span class="rm-lmr-6">1        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-2-4">  <span class="rm-lmr-6">1                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-2-5">  <span class="rm-lmr-6">1          </span></td>
</tr><tr id="TBL-7-3-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-3-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_01.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-3-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-3-3">  <span class="rm-lmr-6">2        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-3-4">  <span class="rm-lmr-6">1                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-3-5">  <span class="rm-lmr-6">2          </span></td>
</tr><tr id="TBL-7-4-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-4-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_01.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-4-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-4-3">  <span class="rm-lmr-6">3        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-4-4">  <span class="rm-lmr-6">1                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-4-5">  <span class="rm-lmr-6">3          </span></td>
</tr><tr id="TBL-7-5-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-5-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_01.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-5-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-5-3">  <span class="rm-lmr-6">4        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-5-4">  <span class="rm-lmr-6">1                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-5-5">  <span class="rm-lmr-6">4          </span></td>
</tr><tr id="TBL-7-6-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-6-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_01.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-6-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-6-3">  <span class="rm-lmr-6">5        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-6-4">  <span class="rm-lmr-6">1                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-6-5">  <span class="rm-lmr-6">5          </span></td>
</tr><tr id="TBL-7-7-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-7-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_01.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-7-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-7-3">  <span class="rm-lmr-6">6        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-7-4">  <span class="rm-lmr-6">1                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-7-5">  <span class="rm-lmr-6">6          </span></td>
</tr><tr id="TBL-7-8-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-8-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_01.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-8-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-8-3">  <span class="rm-lmr-6">7        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-8-4">  <span class="rm-lmr-6">1                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-8-5">  <span class="rm-lmr-6">7          </span></td>
</tr><tr id="TBL-7-9-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-9-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_01.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-9-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-9-3">  <span class="rm-lmr-6">8        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-9-4">  <span class="rm-lmr-6">1                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-9-5">  <span class="rm-lmr-6">8          </span></td>
</tr><tr id="TBL-7-10-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-10-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_01.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-10-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-10-3">  <span class="rm-lmr-6">9        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-10-4">  <span class="rm-lmr-6">1                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-10-5">  <span class="rm-lmr-6">9          </span></td>
</tr><tr id="TBL-7-11-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-11-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_01.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-11-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-11-3">  <span class="rm-lmr-6">10       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-11-4">  <span class="rm-lmr-6">1                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-11-5">  <span class="rm-lmr-6">10         </span></td>
</tr><tr id="TBL-7-12-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-12-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_02.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-12-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-12-3">  <span class="rm-lmr-6">1        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-12-4">  <span class="rm-lmr-6">2                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-12-5">  <span class="rm-lmr-6">11         </span></td>
</tr><tr id="TBL-7-13-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-13-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_02.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-13-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-13-3">  <span class="rm-lmr-6">2        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-13-4">  <span class="rm-lmr-6">2                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-13-5">  <span class="rm-lmr-6">12         </span></td>
</tr><tr id="TBL-7-14-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-14-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_02.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-14-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-14-3">  <span class="rm-lmr-6">3        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-14-4">  <span class="rm-lmr-6">2                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-14-5">  <span class="rm-lmr-6">13         </span></td>
</tr><tr id="TBL-7-15-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-15-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_02.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-15-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-15-3">  <span class="rm-lmr-6">4        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-15-4">  <span class="rm-lmr-6">2                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-15-5">  <span class="rm-lmr-6">14         </span></td>
</tr><tr id="TBL-7-16-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-16-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_02.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-16-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-16-3">  <span class="rm-lmr-6">5        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-16-4">  <span class="rm-lmr-6">2                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-16-5">  <span class="rm-lmr-6">15         </span></td>
</tr><tr id="TBL-7-17-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-17-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_02.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-17-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-17-3">  <span class="rm-lmr-6">6        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-17-4">  <span class="rm-lmr-6">2                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-17-5">  <span class="rm-lmr-6">16         </span></td>
</tr><tr id="TBL-7-18-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-18-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_02.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-18-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-18-3">  <span class="rm-lmr-6">7        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-18-4">  <span class="rm-lmr-6">2                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-18-5">  <span class="rm-lmr-6">17         </span></td>
</tr><tr id="TBL-7-19-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-19-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_02.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-19-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-19-3">  <span class="rm-lmr-6">8        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-19-4">  <span class="rm-lmr-6">2                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-19-5">  <span class="rm-lmr-6">18         </span></td>
</tr><tr id="TBL-7-20-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-20-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_02.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-20-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-20-3">  <span class="rm-lmr-6">9        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-20-4">  <span class="rm-lmr-6">2                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-20-5">  <span class="rm-lmr-6">19         </span></td>
</tr><tr id="TBL-7-21-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-21-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_02.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-21-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-21-3">  <span class="rm-lmr-6">10       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-21-4">  <span class="rm-lmr-6">2                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-21-5">  <span class="rm-lmr-6">20         </span></td>
</tr><tr id="TBL-7-22-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-22-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_03.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-22-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-22-3">  <span class="rm-lmr-6">1        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-22-4">  <span class="rm-lmr-6">3                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-22-5">  <span class="rm-lmr-6">21         </span></td>
</tr><tr id="TBL-7-23-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-23-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_03.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-23-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-23-3">  <span class="rm-lmr-6">2        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-23-4">  <span class="rm-lmr-6">3                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-23-5">  <span class="rm-lmr-6">22         </span></td>
</tr><tr id="TBL-7-24-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-24-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_03.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-24-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-24-3">  <span class="rm-lmr-6">3        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-24-4">  <span class="rm-lmr-6">3                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-24-5">  <span class="rm-lmr-6">23         </span></td>
</tr><tr id="TBL-7-25-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-25-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_03.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-25-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-25-3">  <span class="rm-lmr-6">4        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-25-4">  <span class="rm-lmr-6">3                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-25-5">  <span class="rm-lmr-6">24         </span></td>
</tr><tr id="TBL-7-26-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-26-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_03.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-26-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-26-3">  <span class="rm-lmr-6">5        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-26-4">  <span class="rm-lmr-6">3                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-26-5">  <span class="rm-lmr-6">25         </span></td>
</tr><tr id="TBL-7-27-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-27-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_03.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-27-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-27-3">  <span class="rm-lmr-6">6        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-27-4">  <span class="rm-lmr-6">3                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-27-5">  <span class="rm-lmr-6">26         </span></td>
</tr><tr id="TBL-7-28-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-28-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_03.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-28-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-28-3">  <span class="rm-lmr-6">7        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-28-4">  <span class="rm-lmr-6">3                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-28-5">  <span class="rm-lmr-6">27         </span></td>
</tr><tr id="TBL-7-29-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-29-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_03.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-29-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-29-3">  <span class="rm-lmr-6">8        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-29-4">  <span class="rm-lmr-6">3                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-29-5">  <span class="rm-lmr-6">28         </span></td>
</tr><tr id="TBL-7-30-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-30-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_03.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-30-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-30-3">  <span class="rm-lmr-6">9        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-30-4">  <span class="rm-lmr-6">3                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-30-5">  <span class="rm-lmr-6">29         </span></td>
</tr><tr id="TBL-7-31-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-31-1">  <span class="rm-lmr-6">161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_03.mzML  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-31-2">  <span class="rm-lmr-6">1            </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-31-3">  <span class="rm-lmr-6">10       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-31-4">  <span class="rm-lmr-6">3                     </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-7-31-5">  <span class="rm-lmr-6">30         </span></td></tr></tbody></table>
  </figure>
  </div>


   <div class="table" id="table-5">


   <!-- l. 1777 --><p class="indent">   </p><figure class="float">


<a id="x1-50005r5"></a>
<a id="x1-50006"></a>
   <figcaption class="caption"><span class="id">Table&nbsp;5: Experimental Design 2 </span></figcaption><!-- tex4ht:label?: x1-50005r5  -->
    <table class="tabular" id="TBL-8"><colgroup id="TBL-8-1g"><col id="TBL-8-1"><col id="TBL-8-2"><col id="TBL-8-3"><col id="TBL-8-4"><col id="TBL-8-5"></colgroup><tbody><tr id="TBL-8-1-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-1-1">  <span class="rm-lmr-6">Sample  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-1-2">  <span class="rm-lmr-6">MSstats_Condition  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-1-3">  <span class="rm-lmr-6">MSstats_BioReplicate  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-1-4">  <span class="rm-lmr-6">MSstats_Mixture  </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-1-5">  <span class="rm-lmr-6">LabelName  </span></td>
   </tr><tr id="TBL-8-2-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-2-1">  <span class="rm-lmr-6">1           </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-2-2">  <span class="rm-lmr-6">Norm                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-2-3">  <span class="rm-lmr-6">Norm                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-2-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-2-5">  <span class="rm-lmr-6">126            </span></td>
   </tr><tr id="TBL-8-3-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-3-1">  <span class="rm-lmr-6">2           </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-3-2">  <span class="rm-lmr-6">0.667                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-3-3">  <span class="rm-lmr-6">0.667                        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-3-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-3-5">  <span class="rm-lmr-6">127N          </span></td>
   </tr><tr id="TBL-8-4-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-4-1">  <span class="rm-lmr-6">3           </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-4-2">  <span class="rm-lmr-6">0.125                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-4-3">  <span class="rm-lmr-6">0.125                        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-4-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-4-5">  <span class="rm-lmr-6">127C          </span></td>
   </tr><tr id="TBL-8-5-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-5-1">  <span class="rm-lmr-6">4           </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-5-2">  <span class="rm-lmr-6">0.5                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-5-3">  <span class="rm-lmr-6">0.5                           </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-5-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-5-5">  <span class="rm-lmr-6">128N          </span></td>
   </tr><tr id="TBL-8-6-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-6-1">  <span class="rm-lmr-6">5           </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-6-2">  <span class="rm-lmr-6">1                          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-6-3">  <span class="rm-lmr-6">1                             </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-6-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-6-5">  <span class="rm-lmr-6">128C          </span></td>
   </tr><tr id="TBL-8-7-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-7-1">  <span class="rm-lmr-6">6           </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-7-2">  <span class="rm-lmr-6">0.125                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-7-3">  <span class="rm-lmr-6">0.125                        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-7-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-7-5">  <span class="rm-lmr-6">129N          </span></td>
   </tr><tr id="TBL-8-8-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-8-1">  <span class="rm-lmr-6">7           </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-8-2">  <span class="rm-lmr-6">0.5                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-8-3">  <span class="rm-lmr-6">0.5                           </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-8-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-8-5">  <span class="rm-lmr-6">129C          </span></td>
   </tr><tr id="TBL-8-9-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-9-1">  <span class="rm-lmr-6">8           </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-9-2">  <span class="rm-lmr-6">1                          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-9-3">  <span class="rm-lmr-6">1                             </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-9-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-9-5">  <span class="rm-lmr-6">130N          </span></td>
   </tr><tr id="TBL-8-10-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-10-1">  <span class="rm-lmr-6">9           </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-10-2">  <span class="rm-lmr-6">0.667                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-10-3">  <span class="rm-lmr-6">0.667                        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-10-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-10-5">  <span class="rm-lmr-6">130C          </span></td>
   </tr><tr id="TBL-8-11-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-11-1">  <span class="rm-lmr-6">10          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-11-2">  <span class="rm-lmr-6">Norm                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-11-3">  <span class="rm-lmr-6">Norm                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-11-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-11-5">  <span class="rm-lmr-6">131            </span></td>
   </tr><tr id="TBL-8-12-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-12-1">  <span class="rm-lmr-6">11          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-12-2">  <span class="rm-lmr-6">Norm                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-12-3">  <span class="rm-lmr-6">Norm                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-12-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-12-5">  <span class="rm-lmr-6">126            </span></td>
   </tr><tr id="TBL-8-13-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-13-1">  <span class="rm-lmr-6">12          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-13-2">  <span class="rm-lmr-6">0.667                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-13-3">  <span class="rm-lmr-6">0.667                        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-13-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-13-5">  <span class="rm-lmr-6">127N          </span></td>
   </tr><tr id="TBL-8-14-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-14-1">  <span class="rm-lmr-6">13          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-14-2">  <span class="rm-lmr-6">0.125                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-14-3">  <span class="rm-lmr-6">0.125                        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-14-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-14-5">  <span class="rm-lmr-6">127C          </span></td>
   </tr><tr id="TBL-8-15-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-15-1">  <span class="rm-lmr-6">14          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-15-2">  <span class="rm-lmr-6">0.5                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-15-3">  <span class="rm-lmr-6">0.5                           </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-15-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-15-5">  <span class="rm-lmr-6">128N          </span></td>
   </tr><tr id="TBL-8-16-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-16-1">  <span class="rm-lmr-6">15          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-16-2">  <span class="rm-lmr-6">1                          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-16-3">  <span class="rm-lmr-6">1                             </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-16-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-16-5">  <span class="rm-lmr-6">128C          </span></td>
   </tr><tr id="TBL-8-17-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-17-1">  <span class="rm-lmr-6">16          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-17-2">  <span class="rm-lmr-6">0.125                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-17-3">  <span class="rm-lmr-6">0.125                        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-17-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-17-5">  <span class="rm-lmr-6">129N          </span></td>
   </tr><tr id="TBL-8-18-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-18-1">  <span class="rm-lmr-6">17          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-18-2">  <span class="rm-lmr-6">0.5                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-18-3">  <span class="rm-lmr-6">0.5                           </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-18-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-18-5">  <span class="rm-lmr-6">129C          </span></td>
   </tr><tr id="TBL-8-19-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-19-1">  <span class="rm-lmr-6">18          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-19-2">  <span class="rm-lmr-6">1                          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-19-3">  <span class="rm-lmr-6">1                             </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-19-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-19-5">  <span class="rm-lmr-6">130N          </span></td>
   </tr><tr id="TBL-8-20-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-20-1">  <span class="rm-lmr-6">19          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-20-2">  <span class="rm-lmr-6">0.667                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-20-3">  <span class="rm-lmr-6">0.667                        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-20-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-20-5">  <span class="rm-lmr-6">130C          </span></td>
   </tr><tr id="TBL-8-21-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-21-1">  <span class="rm-lmr-6">20          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-21-2">  <span class="rm-lmr-6">Norm                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-21-3">  <span class="rm-lmr-6">Norm                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-21-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-21-5">  <span class="rm-lmr-6">131            </span></td>
   </tr><tr id="TBL-8-22-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-22-1">  <span class="rm-lmr-6">21          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-22-2">  <span class="rm-lmr-6">Norm                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-22-3">  <span class="rm-lmr-6">Norm                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-22-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-22-5">  <span class="rm-lmr-6">126            </span></td>
   </tr><tr id="TBL-8-23-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-23-1">  <span class="rm-lmr-6">22          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-23-2">  <span class="rm-lmr-6">0.667                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-23-3">  <span class="rm-lmr-6">0.667                        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-23-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-23-5">  <span class="rm-lmr-6">127N          </span></td>
   </tr><tr id="TBL-8-24-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-24-1">  <span class="rm-lmr-6">23          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-24-2">  <span class="rm-lmr-6">0.125                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-24-3">  <span class="rm-lmr-6">0.125                        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-24-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-24-5">  <span class="rm-lmr-6">127C          </span></td>
   </tr><tr id="TBL-8-25-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-25-1">  <span class="rm-lmr-6">24          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-25-2">  <span class="rm-lmr-6">0.5                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-25-3">  <span class="rm-lmr-6">0.5                           </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-25-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-25-5">  <span class="rm-lmr-6">128N          </span></td>
   </tr><tr id="TBL-8-26-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-26-1">  <span class="rm-lmr-6">25          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-26-2">  <span class="rm-lmr-6">1                          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-26-3">  <span class="rm-lmr-6">1                             </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-26-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-26-5">  <span class="rm-lmr-6">128C          </span></td>
   </tr><tr id="TBL-8-27-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-27-1">  <span class="rm-lmr-6">26          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-27-2">  <span class="rm-lmr-6">0.125                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-27-3">  <span class="rm-lmr-6">0.125                        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-27-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-27-5">  <span class="rm-lmr-6">129N          </span></td>
   </tr><tr id="TBL-8-28-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-28-1">  <span class="rm-lmr-6">27          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-28-2">  <span class="rm-lmr-6">0.5                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-28-3">  <span class="rm-lmr-6">0.5                           </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-28-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-28-5">  <span class="rm-lmr-6">129C          </span></td>
   </tr><tr id="TBL-8-29-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-29-1">  <span class="rm-lmr-6">28          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-29-2">  <span class="rm-lmr-6">1                          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-29-3">  <span class="rm-lmr-6">1                             </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-29-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-29-5">  <span class="rm-lmr-6">130N          </span></td>
   </tr><tr id="TBL-8-30-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-30-1">  <span class="rm-lmr-6">29          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-30-2">  <span class="rm-lmr-6">0.667                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-30-3">  <span class="rm-lmr-6">0.667                        </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-30-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-30-5">  <span class="rm-lmr-6">130C          </span></td>
   </tr><tr id="TBL-8-31-"><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-31-1">  <span class="rm-lmr-6">30          </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-31-2">  <span class="rm-lmr-6">Norm                    </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-31-3">  <span class="rm-lmr-6">Norm                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-31-4">  <span class="rm-lmr-6">1                       </span></td><td class="td11" style="white-space:nowrap; text-align:left;" id="TBL-8-31-5">  <span class="rm-lmr-6">131            </span></td></tr></tbody></table>
  </figure>
  </div>

After running the worklfow, the `MSStatsConverter` will convert the OpenMS output in addition with the experimental design to a file (.csv) which can be processed by using `MSstatsTMT`.

#### MSstatsTMT analysis

Here, we depict the analysis by `MSstatsTMT` using a segment of the isobaric analysis workflow (<a href="#figure-25">Fig. 25</a>). The segment is available as {path}`Workflows,MSstatsTMT.knwf`.

(Figure_25)=
|![MSstatsTMT workflow segment](/_images/openms-user-tutorial/isobaric/isobaric_msstatstmt_wf.png)|
|:--:|
|Figure 25: MSstatsTMT workflow segment|

There are two input nodes, the first one takes the result (.csv) from the `MSStatsConverter` and the second a path to the directory where the plots generated by `MSstatsTMT` should be saved. The **R source** node loads the required packages, such as `dplyr` for data wrangling, `MSstatsTMT` for analysis and `MSstats` for plotting. The inputs are further processed in the **R View** node.

Here, the data of the `File Importer` is loaded into **R** using the flow variable [”URI-0”]:

```r
file <- substr(knime.flow.in[["URI-0"]], 6, nchar(knime.flow.in[["URI-0"]]))

MSstatsConverter_OpenMS_out <- read.csv(file)
data <- MSstatsConverter_OpenMS_out
```

The `OpenMStoMSstatsTMTFormat` function preprocesses the OpenMS report and converts it into the required input format for `MSstatsTMT`, by filtering based on unique peptides and measurements in each MS run.

```r
processed.data <- OpenMStoMSstatsTMTFormat(data)
```

Afterwards different normalization steps are performed (global, protein, runs) as well as data imputation by using the msstats method. In addition peptide level data is summarized to protein level data.

```r
quant.data <- proteinSummarization(processed.data,   
                                  method="msstats",

                                  global_norm=TRUE,  
                                  reference_norm=TRUE,

                                  MBimpute = TRUE,  
                                  maxQuantileforCensored = NULL,

                                  remove_norm_channel = TRUE,
                                  remove_empty_channel =  TRUE)
```

There a lot of different possibilities to configure this method please have a look at the MSstatsTMT package for [additional detailed information](http://bioconductor.org/packages/release/bioc/html/MSstatsTMT.html).

The next step is the comparions of the different conditions, here either a pairwise comparision can be performed or a confusion matrix can be created. The goal is to detect and compare the UPS peptides spiked in at different concentrations.

```r
# prepare contrast matrix   
unique(quant.data$Condition)  

comparison<-matrix(c(-1,0,0,1,

                     0,-1,0,1,  
                     0,0,-1,1,

                     0,1,-1,0,  
                     1,-1,0,0), nrow=5, byrow = T)  


# Set the names of each row  
row.names(comparison)<- contrasts <- c("1-0125",

                                       "1-05",  
                                       "1-0667",

                                       "05-0667",  
                                       "0125-05")

# Set the column names
colnames(comparison)<- c("0.125", "0.5", "0.667", "1")
```

The constructed confusion matrix is used in the `groupComparisonTMT` function to test for significant changes in protein abundance across conditions based on a family of linear mixed-effects models in TMT experiments.

```r
data.res <- groupComparisonTMT(data = quant.data,   
                               contrast.matrix = comparison,

                               moderated = TRUE, # do moderated t test

                               adj.method = "BH") # multiple comparison adjustment
data.res <- data.res %>% filter(!is.na(Protein))
```

In the next step the comparison can be plotted using the `groupComparisonPlots` function by `MSstats`.

```r
library(MSstats)  
groupComparisonPlots(data=data.res.mod, type="VolcanoPlot", address=F, which.Comparison = "0125-05", sig = 0.05)
```

Here, we have a example output of the **R View**, which depicts the significant regulated UPS proteins in the comparison of 125 to 05 (<a href="#figure-26">Fig. 26</a>).

(Figure_26)=
|![Volcanoplot of the group comparison regarding 0125 to 05](/_images/openms-user-tutorial/isobaric/isobaric_img_output_knime.png)|
|:--:|
|Figure 26: Volcanoplot of the group comparison regarding 0125 to 05|

All plots are saved to the in the beginning specified output directory in addition.

#### Note

The isobaric analysis does not always has to be performed on protein level, for example for phosphoproteomics studies one is usually interested on the peptide level - in addition inference on peptides with post-translational modification is not straight forward. Here, we present and additonal workflow on peptide level, which can potentially be adapted and used for such cases. Please see {path}`Workflows,IdentificationquantificationisobaricMSstatsTMT`

## References

[^1]: A. Chawade, M. Sandin, J. Teleman, J. Malmström, and F. Levander, Data Processing Has Major Impact on the Outcome of Quantitative Label-Free LC-MS Analysis, Journal of Proteome Research 14(2), 676–687 (2015), PMID: 25407311, arXiv:http://dx.doi.org/10.1021/pr500665j, doi:10.1021/pr500665j. 30

[^2]: M. Choi, Z. F. Eren-Dogu, C. Colangelo, J. Cottrell, M. R. Hoopmann, E. A. Kapp,
S. Kim, H. Lam, T. A. Neubert, M. Palmblad, B. S. Phinney, S. T. Weintraub, B. MacLean, and O. Vitek, ABRF Proteome Informatics Research Group (iPRG)
2015 Study: Detection of Differentially Abundant Proteins in Label-Free Quantitative LC-MS/MS Experiments, J. Proteome Res. 16(2), 945–957 (2017), <a href="https://pubs.acs.org/doi/10.1021/acs.jproteome.6b00881">doi: 10.1021/acs.jproteome.6b00881</a>. 40

[^3]: T. Huang, M. Choi, S. Hao, and O. Vitek, MSstatsTMT: Protein Significance Analysis in shotgun mass spectrometry-based proteomic experiments with tandem
mass tag (TMT) labeling., (2020), <a href="https://bioconductor.org/packages/release/bioc/html/MSstatsTMT.html">doi:10.18129/B9.bioc.MSstatsTMT</a>. 55
