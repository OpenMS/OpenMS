// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:	Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <vector>
#include <map>
#include <set>

namespace OpenMS
{
  class ConsensusMap;
  class FeatureMap;

  /**

  @brief Representation of an experimental design in OpenMS. Instances can be loaded with
         the ExperimentalDesignFile class.

  An experimental design answers the two questions no mzML file can answer by itself:
  <b>which measured value belongs to which piece of biological material</b>, and <b>which of
  those pieces are meant to be compared with each other</b>. Every quantitative OpenMS workflow
  needs both answers, which is why @ref TOPP_ProteomicsLFQ, @ref TOPP_IsobaricWorkflow,
  @ref TOPP_ProteinQuantifier and @ref TOPP_MSstatsConverter all accept a design file.

  The design is a plain TAB-separated text file (conventionally @c design.tsv): a few numeric
  columns that describe the <i>measurement layout</i>, plus any number of free-form columns
  that describe the <i>biology</i>.

  @section ExperimentalDesign_Model The one idea to start from

  <b>One row of the design describes exactly one measured quantity: one channel of one MS file.</b>

  A label-free run carries one channel, so it contributes one row. A TMT10-plex run carries ten
  channels, so it contributes ten rows -- the same file name ten times, with @c Label 1 to 10.
  Every rule below follows from that.

  Each row answers five questions:

  | Column              | Question it answers                                                |
  |---------------------|--------------------------------------------------------------------|
  | @c Spectra_Filepath | Which file was measured?                                           |
  | @c Label            | Which channel inside that file?                                    |
  | @c Fraction         | Which slice of the injected material does this file contain?       |
  | @c Fraction_Group   | Which rows must be combined to obtain one complete measurement?    |
  | @c Sample           | Which piece of biological material was measured here?              |

  These five columns make up the <b>MS file section</b>. The first four are the mechanical part
  of the design and OpenMS validates them strictly (see @ref ExperimentalDesign_Rules).
  @c Sample is the bridge to the <b>sample section</b>: a free-form table of columns -- OpenMS calls them
  <i>factors</i> -- such as condition, biological replicate, genotype or time point. OpenMS never
  interprets the <i>values</i> of those columns; it only compares them to decide which samples
  belong together.

  @section ExperimentalDesign_Glossary Glossary

  These terms are the whole vocabulary. They are easy to mix up because everyday lab language
  uses several of them interchangeably, so each entry says what OpenMS means by it, how it is
  written down, and what OpenMS does with it.

  <b>MS file (run)</b> -- one raw/mzML file, i.e. one LC-MS measurement. Named by
  @c Spectra_Filepath. A file appears in as many rows as it has channels.

  <b>Sample</b> -- the piece of biological material whose abundance you want to know. Technically
  it is a <i>named row of the sample section</i>; @c Sample in the file section points at it by
  name. Two file rows with the <i>same</i> @c Sample value describe measurements of the same
  material and may be pooled by downstream statistics; two rows with different values are
  different material even when all their metadata columns happen to agree. In a labeled
  experiment every channel is its own sample, so one TMT10 file describes ten samples.
  getNumberOfSamples() is the number of rows in the sample section.

  @note @c Sample holds a <b>name</b>, i.e. an arbitrary string -- @c sdrf-pipelines writes
        @c "1", @c "2", ..., but @c "BSA1" or @c "patient_7" are equally valid. In the C++ API,
        MSFileSectionEntry::sample_name is that name while MSFileSectionEntry::sample is the
        <b>zero-based row index</b> of the sample in the sample section. The two rarely coincide;
        do not use one where the other is expected.

  <b>Fraction</b> -- when one sample is separated into several portions <i>before</i> LC-MS
  (high-pH reversed phase, SCX, gel bands, ...), each portion is measured in its own file. Those
  files are fractions of one measurement, <i>not</i> repeats of it: OpenMS aggregates them into a
  single abundance. Numbered from 1; use @c 1 in every row for unfractionated data. Use the
  <i>same</i> numbers in every fraction group so that corresponding fractions line up
  (getFractionToMSFilesMapping() groups by this number across groups).

  <b>Fraction group</b> -- the set of rows that together form one complete measurement of one
  injected mixture: all fractions that have to be combined before a single abundance comes out.
  This is OpenMS' unit of quantification. @ref TOPP_ProteinQuantifier reports one value per
  <tt>(Fraction_Group, Label)</tt> pair and calls that pair an <i>assay</i>. Rules: integers,
  starting at 1, consecutive, no gaps. Unfractionated label-free data has one fraction group per
  file. Measuring the same material a second time means a <b>new</b> fraction group carrying the
  <b>same</b> @c Sample -- that is exactly how a technical replicate is written down.

  <b>Label (channel)</b> -- the multiplexing channel within one file. @c 1 for label-free and
  DIA data; @c 1..n for an n-plex, 1-based and in the canonical channel order of the reagent
  (TMT10-plex: 126&rarr;1, 127N&rarr;2, 127C&rarr;3, 128N&rarr;4, 128C&rarr;5, 129N&rarr;6,
  129C&rarr;7, 130N&rarr;8, 130C&rarr;9, 131&rarr;10; SILAC 2-plex: light&rarr;1, heavy&rarr;2;
  iTRAQ4-plex: 114&rarr;1 ... 117&rarr;4). The number is a <i>position</i>, never a reagent name:
  OpenMS resolves it back to a reagent through the quantification method that produced the data,
  and a mismatch silently swaps channels. Every <tt>(file, label)</tt> pair may occur only once
  in a design.

  <b>Technical replicate</b> -- the same material measured again (re-injection, re-run of the
  same digest). Write it as a second fraction group with the same @c Sample value. OpenMS keeps
  the repeats as separate quantities; collapsing them is the job of the statistics package.

  <b>Biological replicate</b> -- an independent biological unit (another animal, patient,
  culture) measured under a condition. Different biological replicates are different samples, and
  their relationship is declared in the sample section, conventionally in a column called
  @c MSstats_BioReplicate.

  Condition and factors are a property of the <i>sample section</i> and are described in
  @ref ExperimentalDesign_SampleSection.

  @section ExperimentalDesign_Formats The two file formats

  The same design can be written in two ways. Both are TAB-separated; cells are whitespace
  trimmed and lines starting with a hash character (a comment) are ignored.

  <b>One-table format</b> -- everything in a single table. The mandatory columns are
  @c Fraction_Group, @c Fraction and @c Spectra_Filepath; @c Label and @c Sample are optional;
  <i>any further column is taken to be sample metadata</i>. This is the simpler and by far the
  more common form.

  <b>Two-table format</b> -- an MS file section and a sample section, separated by one blank
  line. The file section accepts <i>only</i> @c Fraction_Group, @c Fraction,
  @c Spectra_Filepath, @c Label and @c Sample; any other column there is a parse error. The
  sample section must have a @c Sample column and may carry any number of further columns. Use
  this form when many files share a sample and you do not want to repeat its metadata.

  The format is auto-detected: a file is read as two-table as soon as some line contains a
  @c Sample column header but no @c Fraction_Group column header -- that line is the sample
  header. Otherwise the file is read as one-table.

  @section ExperimentalDesign_Examples Worked examples

  @subsection ExperimentalDesign_Ex1 1. Label-free, unfractionated: 2 conditions, 3 biological replicates

  The simplest design there is: six independent biological units, one file each. Every file is
  its own fraction group and its own sample, and the fraction column is 1 throughout because
  nothing was fractionated.

  | Fraction_Group | Fraction | Spectra_Filepath | Label | Sample | MSstats_Condition | MSstats_BioReplicate |
  |----------------|----------|------------------|-------|--------|-------------------|----------------------|
  | 1              | 1        | ctrl_1.mzML      | 1     | 1      | control           | 1                    |
  | 2              | 1        | ctrl_2.mzML      | 1     | 2      | control           | 2                    |
  | 3              | 1        | ctrl_3.mzML      | 1     | 3      | control           | 3                    |
  | 4              | 1        | drug_1.mzML      | 1     | 4      | treated           | 4                    |
  | 5              | 1        | drug_2.mzML      | 1     | 5      | treated           | 5                    |
  | 6              | 1        | drug_3.mzML      | 1     | 6      | treated           | 6                    |

  6 files, 6 fraction groups, 6 samples, 1 fraction, 1 label, 2 conditions.

  @note The biological replicate numbers run 1..6 and are not restarted per condition. Repeating
        a value under two conditions is how a <i>paired</i> (repeated-measures) design is
        declared -- the same individual measured before and after treatment -- and
        @ref TOPP_MSstatsConverter warns when it sees one, because it changes the statistical
        model. Only reuse a value when the pairing is real.

  @subsection ExperimentalDesign_Ex2 2. Adding technical replicates

  Two patients, each digest injected twice. The repeats carry the <b>same</b> @c Sample (and the
  same metadata) but their <b>own</b> fraction group:

  | Fraction_Group | Fraction | Spectra_Filepath | Label | Sample | MSstats_Condition | MSstats_BioReplicate |
  |----------------|----------|------------------|-------|--------|-------------------|----------------------|
  | 1              | 1        | p1_inj1.mzML     | 1     | 1      | control           | 1                    |
  | 2              | 1        | p1_inj2.mzML     | 1     | 1      | control           | 1                    |
  | 3              | 1        | p2_inj1.mzML     | 1     | 2      | treated           | 2                    |
  | 4              | 1        | p2_inj2.mzML     | 1     | 2      | treated           | 2                    |

  4 files, 4 fraction groups, but only 2 samples. Reusing a sample across fraction groups is
  explicitly allowed and does not merge the groups: OpenMS still reports four quantities and
  leaves it to the statistics package to treat them as repeated measurements.

  @subsection ExperimentalDesign_Ex3 3. Fractionated, label-free

  Two samples, each pre-fractionated into three fractions. Same number of files as example 2, but
  a completely different meaning:

  | Fraction_Group | Fraction | Spectra_Filepath | Label | Sample | MSstats_Condition | MSstats_BioReplicate |
  |----------------|----------|------------------|-------|--------|-------------------|----------------------|
  | 1              | 1        | ctrl_F1.mzML     | 1     | 1      | control           | 1                    |
  | 1              | 2        | ctrl_F2.mzML     | 1     | 1      | control           | 1                    |
  | 1              | 3        | ctrl_F3.mzML     | 1     | 1      | control           | 1                    |
  | 2              | 1        | drug_F1.mzML     | 1     | 2      | treated           | 2                    |
  | 2              | 2        | drug_F2.mzML     | 1     | 2      | treated           | 2                    |
  | 2              | 3        | drug_F3.mzML     | 1     | 2      | treated           | 2                    |

  6 files, 2 fraction groups, 2 samples, 3 fractions. Compare with example 2: there, four files
  produced four quantities; here, six files produce two, because the three files of a fraction
  group are aggregated. Fraction numbers repeat across groups on purpose -- @c ctrl_F2 and
  @c drug_F2 are the same slice of the separation and are aligned with each other.

  @subsection ExperimentalDesign_Ex4 4. TMT10-plex, one mixture, unfractionated

  One file, ten channels, ten samples. The file name repeats on every row; only @c Label and
  @c Sample change.

  | Fraction_Group | Fraction | Spectra_Filepath | Label | Sample | MSstats_Condition | MSstats_BioReplicate | MSstats_Mixture |
  |----------------|----------|------------------|-------|--------|-------------------|----------------------|-----------------|
  | 1              | 1        | tmt_mix1.mzML    | 1     | 1      | control           | 1                    | 1               |
  | 1              | 1        | tmt_mix1.mzML    | 2     | 2      | control           | 2                    | 1               |
  | 1              | 1        | tmt_mix1.mzML    | 3     | 3      | control           | 3                    | 1               |
  | 1              | 1        | tmt_mix1.mzML    | 4     | 4      | control           | 4                    | 1               |
  | 1              | 1        | tmt_mix1.mzML    | 5     | 5      | control           | 5                    | 1               |
  | 1              | 1        | tmt_mix1.mzML    | 6     | 6      | treated           | 6                    | 1               |
  | 1              | 1        | tmt_mix1.mzML    | 7     | 7      | treated           | 7                    | 1               |
  | 1              | 1        | tmt_mix1.mzML    | 8     | 8      | treated           | 8                    | 1               |
  | 1              | 1        | tmt_mix1.mzML    | 9     | 9      | treated           | 9                    | 1               |
  | 1              | 1        | tmt_mix1.mzML    | 10    | 10     | treated           | 10                   | 1               |

  1 file, 1 fraction group, 10 labels, 10 samples. @c Label 1 is TMT channel 126 and @c Label 10
  is channel 131 -- see the glossary for the full order.

  @subsection ExperimentalDesign_Ex5 5. TMT10-plex, two mixtures, fractionated

  Two TMT mixtures, each separated into three fractions: 2 &times; 3 = 6 files, 6 &times; 10 = 60
  rows. Only the first and last rows of each block are shown.

  | Fraction_Group | Fraction | Spectra_Filepath  | Label | Sample | MSstats_Condition | MSstats_BioReplicate | MSstats_Mixture |
  |----------------|----------|-------------------|-------|--------|-------------------|----------------------|-----------------|
  | 1              | 1        | mix1_F1.mzML      | 1     | 1      | control           | 1                    | 1               |
  | ...            | ...      | ...               | ...   | ...    | ...               | ...                  | ...             |
  | 1              | 1        | mix1_F1.mzML      | 10    | 10     | treated           | 10                   | 1               |
  | 1              | 2        | mix1_F2.mzML      | 1     | 1      | control           | 1                    | 1               |
  | ...            | ...      | ...               | ...   | ...    | ...               | ...                  | ...             |
  | 1              | 3        | mix1_F3.mzML      | 10    | 10     | treated           | 10                   | 1               |
  | 2              | 1        | mix2_F1.mzML      | 1     | 11     | control           | 11                   | 2               |
  | ...            | ...      | ...               | ...   | ...    | ...               | ...                  | ...             |
  | 2              | 3        | mix2_F3.mzML      | 10    | 20     | treated           | 20                   | 2               |

  6 files, 2 fraction groups, 3 fractions, 10 labels, 20 samples. Note that a sample repeats
  across the fractions of its own group (same material, different slice) but never across the two
  mixtures (different material, different TMT tube).

  @subsection ExperimentalDesign_Ex6 6. SILAC light/heavy

  Metabolic labeling: two channels per file, two files.

  | Fraction_Group | Fraction | Spectra_Filepath | Label | Sample | MSstats_Condition | MSstats_BioReplicate |
  |----------------|----------|------------------|-------|--------|-------------------|----------------------|
  | 1              | 1        | silac_r1.mzML    | 1     | 1      | control           | 1                    |
  | 1              | 1        | silac_r1.mzML    | 2     | 2      | treated           | 2                    |
  | 2              | 1        | silac_r2.mzML    | 1     | 3      | control           | 3                    |
  | 2              | 1        | silac_r2.mzML    | 2     | 4      | treated           | 4                    |

  @c Label 1 is the light channel, @c Label 2 the heavy one. In a 3-plex, medium is 2 and heavy
  is 3.

  @subsection ExperimentalDesign_Ex7 7. The same design in two-table format

  Example 3 again, with the metadata factored out. The blank line is what separates the sections
  and is mandatory.

  MS file section:

  | Fraction_Group | Fraction | Spectra_Filepath | Label | Sample |
  |----------------|----------|------------------|-------|--------|
  | 1              | 1        | ctrl_F1.mzML     | 1     | 1      |
  | 1              | 2        | ctrl_F2.mzML     | 1     | 1      |
  | 1              | 3        | ctrl_F3.mzML     | 1     | 1      |
  | 2              | 1        | drug_F1.mzML     | 1     | 2      |
  | 2              | 2        | drug_F2.mzML     | 1     | 2      |
  | 2              | 3        | drug_F3.mzML     | 1     | 2      |

  (blank line)

  Sample section:

  | Sample | MSstats_Condition | MSstats_BioReplicate |
  |--------|-------------------|----------------------|
  | 1      | control           | 1                    |
  | 2      | treated           | 2                    |

  @section ExperimentalDesign_SampleSection The sample section: factors and conditions

  Every column of the sample section is a <b>factor</b>. Their names are free-form; OpenMS only
  ever compares values, so @c Genotype, @c Timepoint or @c Dose work exactly like the MSstats
  columns below.

  A <b>condition</b> is a unique combination of the values of all factors <i>except</i>
  @c Sample and <i>except</i> every factor whose column name contains @c "replicate" or
  @c "Replicate". Samples that agree on all remaining factors form one condition -- that is what
  getConditionToSampleMapping(), getSampleToConditionMapping() and
  getPathLabelToConditionMapping() return, and it is how @ref TOPP_ProteomicsLFQ decides which
  runs may be merged.

  @warning The replicate rule is purely a <b>naming convention on the column name</b>. A column
           called @c MSstats_BioReplicate or @c Technical_Replicate is excluded from the
           condition automatically; a column called @c Donor or @c Rep is not, and every distinct
           value in it creates its own condition. Name replicate columns accordingly.

  getUniqueSampleRowToSampleMapping() and getSampleToPrefractionationMapping() apply the weaker
  rule: they ignore only @c Sample and keep the replicate columns, so they group samples whose
  metadata rows are completely identical.

  The conventional MSstats columns, understood by @ref TOPP_MSstatsConverter, are:

  - @c MSstats_Condition: the condition of a sample (e.g. @c control, <tt>"1000 mMol"</tt>). Forwarded
    to MSstats, where it is used to formulate test contrasts.
  - @c MSstats_BioReplicate: identifier of the biological unit. A value that recurs under two
    different conditions declares one unit measured in both, i.e. a paired design; give unrelated
    samples distinct values.
  - @c MSstats_Mixture (isobaric labeling only): identifier of the labeled mixture that was
    measured together in one tube. Samples labeled with different TMT reagents of the same mix
    share a mixture identifier; technical replicates of a mixture keep it as well.

  See the
  <a href="https://www.bioconductor.org/packages/release/bioc/vignettes/MSstats/inst/doc/MSstats.html">MSstats
  manual</a> for what MSstats does with them.

  @section ExperimentalDesign_Rules Rules OpenMS enforces

  Most of these are checked when the file is loaded, and a design that breaks one is rejected.
  They are worth knowing before writing a design by hand:

  - Column headers are matched exactly and are case-sensitive.
  - Fraction groups must be integers, must start at 1 and must be consecutive -- no gaps, no
    renumbering per condition.
  - Each <tt>(Fraction_Group, Fraction, Label)</tt> triple may occur only once.
  - Each <tt>(Spectra_Filepath, Label)</tt> pair may occur only once.
  - In a design that uses a single label (i.e. label-free), each
    <tt>(Fraction_Group, Label)</tt> pair must map to exactly one sample. Two samples in one
    label-free fraction group means the labels are missing or the grouping is wrong.
  - One-table format only: @c Sample becomes mandatory as soon as any @c Label is greater than 1
    -- OpenMS cannot guess which channel is which material. Without a @c Sample column, the
    sample name defaults to the @c Fraction_Group value, which makes every fraction group its own
    sample.
  - Two-table format only: every @c Sample value used in the file section must exist in the
    sample section. Omitting the @c Sample column from a two-table file section is possible but
    almost never what you want: OpenMS then looks for samples literally named
    <tt>"Fraction group 1"</tt>, <tt>"Fraction group 2"</tt>, ...
  - A relative @c Spectra_Filepath is resolved first against the directory of the design file,
    then against the current working directory; if neither exists the string is kept as written.
    Tools that need the spectra present (@c require_spectra_files) fail instead.

  @section ExperimentalDesign_SDRF Relation to SDRF-Proteomics

  <a href="https://github.com/bigbio/proteomics-sample-metadata">SDRF-Proteomics</a> is the
  HUPO-PSI sample-and-data-relationship format used by PRIDE and quantms. It is richer than an
  OpenMS design -- it also carries organism, disease, instrument and search settings -- and it is
  the recommended source of truth: <tt>parse_sdrf convert-openms -s sdrf.tsv</tt> from
  <a href="https://github.com/bigbio/sdrf-pipelines">sdrf-pipelines</a> generates the design file
  described here. The conversion maps:

  | OpenMS design         | SDRF-Proteomics                                        | Notes                                                                                     |
  |-----------------------|--------------------------------------------------------|-------------------------------------------------------------------------------------------|
  | @c Spectra_Filepath   | <tt>comment[data file]</tt>                                  | the converter can rewrite the extension, e.g. @c .raw to @c .mzML                          |
  | @c Fraction           | <tt>comment[fraction identifier]</tt>                        | @c 1 when absent or "not available"                                                        |
  | @c Label              | @c comment[label]                                      | the ontology term becomes its 1-based channel index: "label free sample"&rarr;1, TMT126&rarr;1 ... TMT131&rarr;10, "SILAC light"/"SILAC heavy"&rarr;1/2 |
  | @c Fraction_Group     | <tt>source name</tt> + <tt>comment[technical replicate]</tt>  | one group per (biological source, technical replicate) -- which is why re-injections get their own group |
  | @c Sample             | <tt>source name</tt>                                    | a source name ending in "sample N" contributes N, otherwise names are numbered in order of appearance |
  | @c MSstats_Condition  | the <tt>factor value[...]</tt> columns, concatenated    | falls back to the non-redundant <tt>characteristics[...]</tt> columns, then to the source name   |
  | @c MSstats_BioReplicate | derived from <tt>source name</tt>                     | one value per biological source; <i>not</i> copied from <tt>characteristics[biological replicate]</tt> |
  | @c MSstats_Mixture    | derived per labeled file and sample                    | isobaric designs only                                                                      |

  Two SDRF terms translate less directly than their names suggest. SDRF's
  <tt>comment[technical replicate]</tt> is folded into @c Fraction_Group, because OpenMS expresses
  technical replication as "same sample, new fraction group" rather than as a column. And SDRF
  distinguishes <tt>source name</tt> (the material) from <tt>assay name</tt> (the measurement),
  where OpenMS uses @c Sample and the file/label pair for the same distinction.

  @section ExperimentalDesign_QPX Relation to QPX (quantms.io)

  QPX ("Quantitative Proteomics eXchange") is the Parquet exchange format written by
  @ref TOPP_ProteomicsLFQ and @ref TOPP_IsobaricWorkflow via their @c out_qpx option. Its
  sample metadata comes from SDRF, so a design exported to QPX has to line up with the terms
  above:

  | OpenMS design                              | QPX                                              | Notes                                                                                   |
  |--------------------------------------------|--------------------------------------------------|-----------------------------------------------------------------------------------------|
  | basename of @c Spectra_Filepath, no extension | @c run_file_name                               | primary-key component of the psm, feature and pg views; see ArrowIOHelpers::qpxRunFileName() |
  | @c Fraction_Group                          | @c grouped_runs of the pg view                   | QPX has no fraction-group number: the group is represented by listing its run names      |
  | @c Fraction                                | @c run.fraction                                  |                                                                                          |
  | @c Label                                   | <tt>intensities[].label</tt>, joined to <tt>run.samples[].label</tt> | as a canonical channel token (@c "LFQ", @c "TMT126", <tt>"SILAC heavy"</tt>), not as the number; see ArrowIOHelpers::qpxIntensityLabel() |
  | @c Sample                                  | <tt>run.samples[].sample_accession</tt>                | foreign key into @c sample.parquet, which is keyed by the SDRF source name               |
  | replicate factors of the sample section    | <tt>run.samples[].biological_replicate</tt> / <tt>technical_replicate</tt> |                                                                          |

  @note QPX requires a fraction group to be <b>rectangular</b>: every label the group publishes
        must exist in every one of its runs, so that <tt>(any run of the group, label)</tt>
        resolves to one sample. Ragged designs, and designs that put one run into several
        fraction groups, are refused by the exporter.

  @section ExperimentalDesign_Pitfalls Common pitfalls

  - <b>Renumbering fraction groups per condition.</b> They are global: 1, 2, 3, ... across the
    whole design. Restarting at 1 for the second condition is a common load error.
  - <b>Confusing fractions with replicates.</b> Fractions are added together (one quantity out of
    several files); replicates are kept apart (several quantities). Compare examples 2 and 3.
  - <b>Forgetting that a labeled file needs one row per channel.</b> A TMT10 design with one row
    per file quantifies one channel and silently ignores nine.
  - <b>Giving a replicate column a name without "replicate" in it.</b> It then counts as part of
    the condition and splits every group of replicates into separate conditions.
  - <b>Reusing a biological replicate id across conditions by accident.</b> That declares a
    paired design; @ref TOPP_MSstatsConverter warns about it.
  - <b>Using a comma-separated file.</b> The format is TAB-separated, despite the @c .tsv files
    often being edited in a spreadsheet program.

  @see ExperimentalDesignFile for loading, and fromConsensusMap() / fromFeatureMap() /
       fromIdentifications() for deriving a design when no file is available.

  @ingroup Metadata

  **/

  class OPENMS_DLLAPI ExperimentalDesign
  {

  public:
    /**
      @brief One row of the MS file section: one quantitative channel of one MS file

      This is the in-memory form of a row of the design table described on ExperimentalDesign.
      It links a single quantitative value back to the file it came from and supports:
       - multiplexed/labeled data via specification of the quantified label
       - multiple fractions via specification of the:
         - fraction index (e.g., 1..10 if ten fractions were measured)
         - fraction_group to trace which fractions belong together
    */
    class OPENMS_DLLAPI MSFileSectionEntry
    {
    public:
      MSFileSectionEntry() = default;
      /// Fraction group: the rows that are combined into one quantity. 1-based and consecutive
      /// over the whole design. A fraction group is OpenMS' unit of quantification, so a second
      /// measurement of the same material is a new group carrying the same @c sample.
      unsigned fraction_group = 1;
      /// Fraction 1..m, mandatory, 1 if not set (and 1 everywhere for unfractionated data).
      /// The same number in two fraction groups denotes corresponding fractions.
      unsigned fraction = 1;
      std::string path = "UNKNOWN_FILE"; ///< file name, mandatory
      /// The 1-based channel position within @c path (1 for label-free; e.g. 1..10 for
      /// TMT10plex in canonical channel order 126, 127N, ..., 131). A position, not a reagent name.
      unsigned label = 1;
      /// Zero-based ROW INDEX into the sample section; allows grouping by sample.
      /// This is not the @c Sample column value -- that is @c sample_name.
      unsigned sample = 0;
      /// The @c Sample column value, i.e. the sample's NAME (an arbitrary string such as "1" or
      /// "BSA1"). Used to look the sample row up by name rather than by index.
      std::string sample_name = "0";
    };

    /**
      @brief The sample section: one named row per sample, one column per factor

      Holds the free-form metadata half of the design (see @ref ExperimentalDesign_SampleSection).
      Every column is a <i>factor</i>; OpenMS never interprets factor values, it only compares
      them to group samples into conditions.
    */
    class OPENMS_DLLAPI SampleSection
    {
    public:

      SampleSection() = default;

      SampleSection(
        const std::vector< std::vector < std::string > >& content,
        const std::map< std::string, Size >& sample_to_rowindex,
        const std::map< std::string, Size >& columnname_to_columnindex
      );

      /// Get set of all sample NAMES (the @c Sample column values) present in the sample section
      std::set< std::string > getSamples() const;

      /// Add a sample as the last row
      void addSample(const std::string& sample, const std::vector<std::string>& content = {});

      // TODO should it include the Sample ID column or not??
      /// Get set of all factors (column names) that were defined for the sample section.
      /// Note that this currently includes the @c Sample column itself; the mapping functions of
      /// ExperimentalDesign therefore drop it explicitly before comparing rows.
      std::set< std::string > getFactors() const;

      /// Checks whether the sample section has a row for a sample name
      bool hasSample(const std::string& sample) const;

      /// Checks whether Sample Section has a specific factor (i.e. column name)
      bool hasFactor(const std::string &factor) const;

      /// Returns value of factor for given sample NAME and factor name
      std::string getFactorValue(const std::string& sample_name, const std::string &factor) const;

      /// Returns value of factor for given sample ROW INDEX (zero-based) and factor name
      std::string getFactorValue(unsigned sample_idx, const std::string &factor) const;

      /// Returns column index of factor. This is a storage detail, not a property of the design:
      /// a one-table design orders the sample columns alphabetically, a two-table design keeps
      /// the order of the file. Address factors by name unless the layout itself matters.
      Size getFactorColIdx(const std::string &factor) const;

      /// Returns the name/ID of the sample for a zero-based row index. Not the row index itself
      std::string getSampleName(unsigned sample_row) const;

      /// Returns the zero-based row index in the sample section for a sample name/ID
      unsigned getSampleRow(const std::string& sample) const;

      /// returns the number of entries in content_ member
      Size getContentSize() const;

    private:

      // The entries of the Sample Section, filled while parsing
      // the Experimental Design File
      std::vector< std::vector < std::string > > content_;

      // Maps the Sample Entry name to the row where the sample
      // appears in the Sample section, its sample index
      std::map< std::string, Size > sample_to_rowindex_;

      // Maps the column name of the SampleSection to the
      // Index of the column
      std::map< std::string, Size > columnname_to_columnindex_;
    };

    using MSFileSection = std::vector<MSFileSectionEntry>;

    // Experimental Design c'tors
    ExperimentalDesign() = default;

    ExperimentalDesign(const MSFileSection& msfile_section, const SampleSection& sample_section);

    const MSFileSection& getMSFileSection() const;

    void setMSFileSection(const MSFileSection& msfile_section);

    // Returns the Sample Section of the experimental design file
    const ExperimentalDesign::SampleSection& getSampleSection() const;

    void setSampleSection(const SampleSection& sample_section);

    /// returns a map from a sample section row to sample id for clustering
    /// duplicate sample rows (e.g. to find all fractions of the same "sample").
    /// Rows are compared over ALL factors except @c Sample -- replicate columns included -- so
    /// samples are grouped only when their metadata is completely identical. This is the weaker
    /// of the two grouping rules; see getConditionToSampleMapping() for the other one.
    std::map<std::vector<std::string>, std::set<std::string>> getUniqueSampleRowToSampleMapping() const;

    /// uses getUniqueSampleRowToSampleMapping to get the reversed map
    /// mapping sample ID to a real unique sample.
    /// Keyed by sample NAME (the Sample column value), not by the zero-based sample row index.
    std::map<std::string, unsigned> getSampleToPrefractionationMapping() const;

    /// return fraction index to file paths (ordered by fraction_group)
    //TODO this probably needs a basename parameter to be fully compatible with the other mappings!! Implicit full path.
    std::map<unsigned int, std::vector<std::string> > getFractionToMSFilesMapping() const;

    /// return vector of filepath/label combinations that share the same conditions after removing
    /// replicate columns in the sample section (e.g. for merging across replicates)
    //TODO this probably needs a basename parameter to be fully compatible with the other mappings!! Implicit full path.
    std::vector<std::vector<std::pair<std::string, unsigned>>> getConditionToPathLabelVector() const;

    /**
      @brief return a condition to Sample index mapping

      A condition is the unique combination of the sample section's factor values, ignoring the
      @c Sample column and every column whose NAME contains @c "replicate" or @c "Replicate".
      That naming convention is the only thing that marks a column as a replicate column: a
      column called @c MSstats_BioReplicate is excluded, one called @c Donor is not and would
      split the replicates it names into separate conditions.

      @return condition (the retained factor values, in the alphabetical order of their column
              names) to the zero-based sample row indices that share it
    */
    std::map<std::vector<std::string>, std::set<unsigned>> getConditionToSampleMapping() const;

   /*
    *   The (Path, Label) tuples in the experimental design have to be unique, so we can map them
    *   uniquely to the sample number, fraction number, and fraction_group number
    */

    /// return <file_path, label> to prefractionation mapping (a prefractionation group is a unique combination of
    /// all columns in the sample section, except for replicates.
    std::map< std::pair< std::string, unsigned >, unsigned> getPathLabelToPrefractionationMapping(bool use_basename_only) const;

    /// return <file_path, label> to condition mapping (a condition is a unique combination of all columns in the
    /// sample section, except for replicates.
    std::map< std::pair< std::string, unsigned >, unsigned> getPathLabelToConditionMapping(bool use_basename_only) const;

    /// return Sample name to condition mapping (a condition is a unique combination of all columns in the
    /// sample section, except for replicates. Numbering of conditions is alphabetical due to map.
    /// Keyed by sample NAME, like getSampleToPrefractionationMapping(); it previously used the
    /// stringified zero-based sample row index, which only matched designs OpenMS inferred itself.
    std::map<std::string, unsigned> getSampleToConditionMapping() const;

    /// return <file_path, label> to sample index mapping
    std::map< std::pair< std::string, unsigned >, unsigned> getPathLabelToSampleMapping(bool use_basename_only) const;

    /// return <file_path, label> to fraction mapping
    std::map< std::pair< std::string, unsigned >, unsigned> getPathLabelToFractionMapping(bool use_basename_only) const;

    /// return <file_path, label> to fraction_group mapping
    std::map< std::pair< std::string, unsigned >, unsigned> getPathLabelToFractionGroupMapping(bool use_basename_only) const;

    /// @return the number of samples measured (= number of rows in the sample section).
    /// NOT the highest sample index: MSFileSectionEntry::sample is zero-based, so the highest
    /// index is one less than this.
    unsigned getNumberOfSamples() const;

    /// @return the number of distinct fraction indices used anywhere in the design.
    /// Meaningful only when all fraction groups use the same fraction numbering; see
    /// sameNrOfMSFilesPerFraction() to check that.
    unsigned getNumberOfFractions() const;

    /// @return the number of labels per file (= the highest label index, i.e. the plex size;
    /// 1 for a label-free design)
    unsigned getNumberOfLabels() const;

    /// @return the number of distinct MS file paths in the design. A multiplexed file
    /// contributes one file but several rows.
    unsigned getNumberOfMSFiles() const;

    /// @return the number of distinct fraction_groups.
    /// Allows to group fraction ids and source files
    unsigned getNumberOfFractionGroups() const;

    /// @return the zero-based sample row index quantified in @p fraction_group at @p label
    /// @throw Exception::ElementNotFound if the design has no such combination
    unsigned getSample(unsigned fraction_group, unsigned label = 1);

    /// @return whether we have a fractionated design
    /// This is the case if more than one distinct fraction index is used in the design
    bool isFractionated() const;

    /// filters the MSFileSection to only include a given subset of files whose basenames
    /// are given with @p bns
    /// @return number of files that have been filtered
    Size filterByBasenames(const std::set<std::string>& bns);

    /// @returns whether all fraction groups have the same number of fractions
    bool sameNrOfMSFilesPerFraction() const;

    /// Extract experimental design from consensus map
    static ExperimentalDesign fromConsensusMap(const ConsensusMap& c);

    /**
      @brief Write this design's fraction structure onto a ConsensusMap's column headers

      The inverse of fromConsensusMap(): stamps @c fraction_group, @c fraction and
      @c sample_name as meta values on every column header whose @c (basename, label) pair
      matches a design row. Headers are matched on the 1-based label, so this works for both
      label-free maps (one header per file) and multiplexed ones (one header per file and
      channel).

      Without this the fraction structure lives only in the design object, and consumers that
      read the map -- exporters, converters -- cannot tell two fractions of one sample from two
      independent runs.

      @param[in,out] cmap Map whose column headers are annotated in place
      @return Number of column headers with no matching design row, left unannotated. Non-zero
              means the design does not describe the whole map; callers should report it.
    */
    Size annotateColumnHeaders(ConsensusMap& cmap) const;

    /// Extract experimental design from feature map
    static ExperimentalDesign fromFeatureMap(const FeatureMap& f);

    /// Extract experimental design from identifications
    static ExperimentalDesign fromIdentifications(const std::vector<ProteinIdentification>& proteins);
    //TODO create another overload here, that takes two enums outerVec and innerVec with entries Replicate, Fraction, Sample

    private:
    // sample row index -> sample name; the reverse of SampleSection's name->row store
    std::map<Size, std::string> sampleRowToName_() const;

    // MS filename column, optionally trims to basename
    std::vector< std::string > getFileNames_(bool basename) const;

    // returns label column
    std::vector<unsigned> getLabels_() const;

    // returns fraction column
    std::vector<unsigned> getFractions_() const;

    /// Generic Mapper (Path, Label) -> f(row)
    std::map< std::pair< std::string, unsigned >, unsigned> pathLabelMapper_(
        bool,
        unsigned (*f)(const ExperimentalDesign::MSFileSectionEntry&)) const;

    // sort to obtain the default order
    void sort_();

    template<typename T>
    static void errorIfAlreadyExists(std::set<T> &container, T &item, const std::string &message);

    // basic consistency checks
    void isValid_();

    MSFileSection msfile_section_;
    SampleSection sample_section_;
  };
}

