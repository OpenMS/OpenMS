Identification and Quantification
=================================

While the combination of liquid chromatography and mass spectrometry can ease the process of characterising molecules of interest, further techniques are required to easily identify and quantify these molecules. This section discusses both labeled and label-free quantification techniques.

## Labeling

Relative quantification is one strategy where one sample is chemically treated and compared to another sample without treatment. This section discusses a particular relative quanitification technique called **labeling** or **stable isotope labeling** which involves the addition of isotopes to one sample. An isotope of an element behaves the same chemically but has a different mass. Stable isotope labeling is used in mass spectrometry so that scientists can easily identify proteins and metabolites.

Two types of stable isotope labeling exist: chemical labeling and metabolic labeling.

### Chemical labeling

During chemical labeling, the label is attached at specific functional groups in a molecule like the N-terminus of a peptide or specific side chains.

Chemical labeling occurs late in the process, therefore experiments that incorporate this technique are not highly reproducible.

#### Isobaric labeling

Isobaric labeling, is a technique where peptides and proteins are labeled with chemical groups that have an identical mass, but vary in terms of of distribution of heavy isotopes in their structure.

<div class="admonition video">
<p class="admonition-title">**Video**</p>
For more information on isobaric labeling, view the following links:
<ul>
<li><a href="https://timms.uni-tuebingen.de:/tp/UT_20141118_002_cpm_0001?t=1108.15">Video 1</a>
</li>
<li><a href="https://timms.uni-tuebingen.de:/tp/UT_20141202_002_cpm_0001?t=311.78">Video 2</a>
</li>
<ul>
</div>

OpenMS contains tools that analyze data from isobaric labeling experiments.

### Metabolic labeling

During metabolic labeling, the organism is 'fed' with labeled metabolites. Metabolites include but are not limited to amino acids, nitrogen sources and glucose. Unlike chemical labeling, metabolic labeling occurs early in the study. Therefore, experiments that incorporate metabolic labeling are highly reproducible.

#### Stable Isotope Labeling with Amino Aids in Cell Culture (SILAC)

In SILAC, the labeled amino acids are fed to the cell culture. The labels are integrated into the proteins after a period. The labeled sample is then compared with the unlabeled sample.

OpenMS contains tools that analyze data from SILAC experiments.

<div class="admonition video">
<p class="admonition-title">**Video**</p>
For more information on SILAC, view the following links:
<ul>
<li><a href="https://timms.uni-tuebingen.de:/tp/UT_20141118_002_cpm_0001?t=18.25">Video 1</a></li>
<li><a href="https://timms.uni-tuebingen.de:/tp/UT_20141202_001_cpm_0001?t=540.13">Video 2</a></li>
</ul>
</div>

## Label-free quantification (LFQ)

LFQ is a cheap and natural method of quantifying molecules of interest. As the name suggests, no labeling of molecules is involved.

LFQ includes the following steps:

1. **Conduct replicate experiments**.
2. **Generate LC-MS maps** for each experiment.
3. **Find features** in all LC-MS maps. A {term}`feature` is a collection of peaks that belong to a chemical compound.
4. **Align maps** to address shifts in retention times.
5. **Match corresponding features** in different maps. We refer to this as **grouping** or **linking**.
6. **Identify feature groups**, called {term}`consensus features <consensus feature>`.
7. **Quantify consensus features**.

<div class="admonition video">
<p class="admonition-title">**Video**</p>
For more information on LFQ, [view this video](https://timms.uni-tuebingen.de:/tp/UT_20141118_002_cpm_0001?t=2115.00).
For more information on the steps involved in LFQ, [view this video](https://timms.uni-tuebingen.de:/tp/UT_20141118_002_cpm_0001?t=2230.18)
</div>

### Feature finding

Feature finding is method for identifying all peaks belonging to a chemical compound. Feature finding involves the following steps:

1. **Extension** where we collect all data points we think belong to the peptide.
2. **Refinement** where we remove peaks that we think do not belong to the peptide.
3. **Fit an optimal model** to the isolated peaks.

The above steps are iterative; we repeat these steps until no improvement can be made to the model.

OpenMS contains a number of feature finding algorithms.

<div class="admonition video">
<p class="admonition-title">**Video**</p>
For more information on feature finding, [view this video](https://timms.uni-tuebingen.de:/tp/UT_20141118_002_cpm_0001?t=2670.44).
</div>