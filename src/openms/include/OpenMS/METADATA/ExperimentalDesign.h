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

  An experimental design maps quantitative values to the biological material they were measured
  from, and records the metadata needed to compare those measurements. @ref TOPP_ProteomicsLFQ,
  @ref TOPP_IsobaricWorkflow, @ref TOPP_ProteinQuantifier and @ref TOPP_MSstatsConverter all
  take one.

  It is a TAB-separated text file (conventionally @c design.tsv) with a few numeric columns
  describing the measurement layout plus any number of free-form columns describing the biology.

  @section ExperimentalDesign_Model Structure

  <b>One row describes one measured quantity: one channel of one MS file.</b> A label-free run
  has one channel and contributes one row; a TMT10-plex run has ten and contributes ten rows --
  the same file name ten times, with @c Label 1 to 10.

  Per row:

  | Column              | Meaning                                                            |
  |---------------------|--------------------------------------------------------------------|
  | @c Spectra_Filepath | The file that was measured                                         |
  | @c Label            | The channel inside that file                                       |
  | @c Fraction         | Which slice of the injected material the file contains             |
  | @c Fraction_Group   | Which rows are combined into one complete measurement              |
  | @c Sample           | The biological material measured in this row                       |

  These five columns make up the <b>MS file section</b>. The first four are validated strictly
  (see @ref ExperimentalDesign_Rules). @c Sample links to the <b>sample section</b>: a free-form
  table of columns -- OpenMS calls them <i>factors</i> -- such as condition, biological replicate,
  genotype or time point. OpenMS does not interpret factor values; it only compares them to
  decide which samples belong together.

  @section ExperimentalDesign_Glossary Glossary

  Lab usage of these terms varies, so each entry gives the OpenMS meaning, how it is written in
  the file, and what OpenMS does with it.

  <b>MS file (run)</b> -- one raw/mzML file, i.e. one LC-MS measurement. Named by
  @c Spectra_Filepath. A file appears in as many rows as it has channels.

  <b>Sample</b> -- the biological material a quantity is measured from; concretely, a named row
  of the sample section that @c Sample in the file section points at. Two file rows with the same
  @c Sample value describe measurements of the same material; rows with different values are
  different material, even if all their metadata columns agree. Same-sample rows within one
  fraction group are fractions of one measurement and are combined into one quantity; same-sample
  rows in different fraction groups stay separate quantities, never pooled by OpenMS. In a
  labeled experiment each channel normally
  names a different sample, so one TMT10 file usually describes ten, but repeating a @c Sample
  across channels or mixtures is legal and is how a bridge/reference channel is declared.
  getNumberOfSamples() returns the number of rows in the sample section.

  @note @c Sample holds a <b>name</b>, i.e. an arbitrary string -- @c sdrf-pipelines writes
        @c "1", @c "2", ..., but @c "BSA1" or @c "patient_7" are equally valid. In the C++ API,
        MSFileSectionEntry::sample_name is that name while MSFileSectionEntry::sample is the
        <b>zero-based row index</b> of the sample in the sample section. The two rarely coincide;
        do not use one where the other is expected.

  <b>Fraction</b> -- one portion of a sample that was separated before LC-MS (high-pH reversed
  phase, SCX, gel bands, ...) and measured in its own file. Fractions are parts of one
  measurement rather than repeats of it: tools that report one value per assay
  (@ref TOPP_ProteinQuantifier, mzTab, QPX) combine them -- summed by default, or reduced to one
  fraction with @c fractions:aggregate @c best. The MSstats export does not: it emits one row per
  fraction and lets MSstats summarize.
  Numbered from 1; use @c 1 in every row for unfractionated data. Use the same numbers in every
  fraction group so corresponding fractions line up (getFractionToMSFilesMapping() groups by this
  number across groups).

  <b>Fraction group</b> -- the rows that form one complete measurement of one injected mixture:
  all fractions combined before a quantity is reported. The unit a value is actually reported for
  is the pair <tt>(Fraction_Group, Label)</tt> -- @ref TOPP_ProteinQuantifier and mzTab call it an
  <i>assay</i>, the QPX exporters a <i>quantification unit</i>. A label-free fraction group
  therefore yields one quantity, a TMT10 group ten.
  Must be integers starting at 1, consecutive, no gaps.
  Unfractionated label-free data has one fraction group per file. Measuring the same material
  again means a new fraction group with the same @c Sample; see technical replicate below.

  <b>Label (channel)</b> -- the multiplexing channel within one file. @c 1 for label-free and
  DIA data; <tt>1..n</tt> for an n-plex, 1-based in the canonical channel order of the reagent
  (TMT10-plex: 126&rarr;1, 127N&rarr;2, 127C&rarr;3, 128N&rarr;4, 128C&rarr;5, 129N&rarr;6,
  129C&rarr;7, 130N&rarr;8, 130C&rarr;9, 131&rarr;10; SILAC 2-plex: light&rarr;1, heavy&rarr;2;
  iTRAQ4-plex: 114&rarr;1 ... 117&rarr;4). The value is a channel position, not a reagent name:
  OpenMS resolves it to a reagent through the quantification method that produced the data, so a
  mismatch swaps channels without any error. Every <tt>(file, label)</tt> pair may occur only
  once in a design.

  <b>Technical replicate</b> -- the same material measured again (re-injection, re-run of the
  same digest). Written as a second fraction group with the same @c Sample value. OpenMS keeps
  the repeats as separate quantities; collapsing them is left to the statistics package.

  <b>Biological replicate</b> -- an independent biological unit (another animal, patient,
  culture) measured under a condition. Different biological replicates are different samples;
  the relationship is declared in the sample section, conventionally in a column called
  @c MSstats_BioReplicate.

  Conditions and factors belong to the sample section, see @ref ExperimentalDesign_SampleSection.

  @section ExperimentalDesign_Formats The two file formats

  Both are TAB-separated (tabs only -- aligning columns with spaces does not work). While
  parsing, cells are whitespace trimmed and lines starting with a hash character are ignored.

  <b>One-table format</b> -- a single table. Mandatory columns are @c Fraction_Group,
  @c Fraction and @c Spectra_Filepath; @c Label and @c Sample are optional; any further column
  is taken to be sample metadata. This is the more common form.

  <b>Two-table format</b> -- an MS file section and a sample section, separated by at least one
  blank line. The file section accepts only @c Fraction_Group, @c Fraction, @c Spectra_Filepath,
  @c Label and @c Sample; any other column there is a parse error. The sample section must have
  a @c Sample column and may carry any number of further columns. Useful when many files share a
  sample, to avoid repeating its metadata.

  The format is auto-detected, and the detector is cruder than the parsers: it scans <i>every</i>
  line, not just headers, and reads the file as two-table as soon as one line has exactly one cell
  equal to @c Sample and no cell equal to @c Fraction_Group. It does not trim cells or skip
  comment lines, so a sample header whose cells carry padding is missed and a commented-out line
  still counts.

  @section ExperimentalDesign_Examples Worked examples

  @subsection ExperimentalDesign_Ex1 1. Label-free, unfractionated: 2 conditions, 3 biological replicates

  Six independent biological units, one file each. Every file is its own fraction group and its
  own sample; the fraction column is 1 throughout because nothing was fractionated.

  | Fraction_Group | Fraction | Spectra_Filepath | Label | Sample | MSstats_Condition | MSstats_BioReplicate |
  |----------------|----------|------------------|-------|--------|-------------------|----------------------|
  | 1              | 1        | ctrl_1.mzML      | 1     | 1      | control           | 1                    |
  | 2              | 1        | ctrl_2.mzML      | 1     | 2      | control           | 2                    |
  | 3              | 1        | ctrl_3.mzML      | 1     | 3      | control           | 3                    |
  | 4              | 1        | drug_1.mzML      | 1     | 4      | treated           | 4                    |
  | 5              | 1        | drug_2.mzML      | 1     | 5      | treated           | 5                    |
  | 6              | 1        | drug_3.mzML      | 1     | 6      | treated           | 6                    |

  6 files, 6 fraction groups, 6 samples, 1 fraction, 1 label, 2 conditions.

  @note The biological replicate numbers run 1..6 and are not restarted per condition. A value
        repeated under two conditions declares a <i>paired</i> (repeated-measures) design, e.g.
        the same individual before and after treatment. It changes the statistical model, so
        @ref TOPP_MSstatsConverter warns about it; only reuse a value when the pairing is real.

  @subsection ExperimentalDesign_Ex2 2. Adding technical replicates

  Two patients, each digest injected twice. The repeats carry the same @c Sample and metadata,
  but their own fraction group:

  | Fraction_Group | Fraction | Spectra_Filepath | Label | Sample | MSstats_Condition | MSstats_BioReplicate |
  |----------------|----------|------------------|-------|--------|-------------------|----------------------|
  | 1              | 1        | p1_inj1.mzML     | 1     | 1      | control           | 1                    |
  | 2              | 1        | p1_inj2.mzML     | 1     | 1      | control           | 1                    |
  | 3              | 1        | p2_inj1.mzML     | 1     | 2      | treated           | 2                    |
  | 4              | 1        | p2_inj2.mzML     | 1     | 2      | treated           | 2                    |

  4 files, 4 fraction groups, 2 samples. Reusing a sample across fraction groups is allowed and
  does not merge them: OpenMS reports four quantities and leaves it to the statistics package to
  treat them as repeated measurements.

  @subsection ExperimentalDesign_Ex3 3. Fractionated, label-free

  Two samples again, as in example 2, but six files instead of four:

  | Fraction_Group | Fraction | Spectra_Filepath | Label | Sample | MSstats_Condition | MSstats_BioReplicate |
  |----------------|----------|------------------|-------|--------|-------------------|----------------------|
  | 1              | 1        | ctrl_F1.mzML     | 1     | 1      | control           | 1                    |
  | 1              | 2        | ctrl_F2.mzML     | 1     | 1      | control           | 1                    |
  | 1              | 3        | ctrl_F3.mzML     | 1     | 1      | control           | 1                    |
  | 2              | 1        | drug_F1.mzML     | 1     | 2      | treated           | 2                    |
  | 2              | 2        | drug_F2.mzML     | 1     | 2      | treated           | 2                    |
  | 2              | 3        | drug_F3.mzML     | 1     | 2      | treated           | 2                    |

  6 files, 2 fraction groups, 2 samples, 3 fractions. In example 2 four files produced four
  quantities; here six files produce two, because the files of a fraction group are aggregated.
  Fraction numbers repeat across groups deliberately: @c ctrl_F2 and @c drug_F2 are the same
  slice of the separation and are aligned with each other.

  @subsection ExperimentalDesign_Ex4 4. TMT10-plex, one mixture, unfractionated

  One file, ten channels, ten samples. The file name repeats on every row; of the measurement
  columns only @c Label changes, and each channel names its own sample.

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
  rows, of which a few representative ones are shown.

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

  6 files, 2 fraction groups, 3 fractions, 10 labels, 20 samples. A sample repeats across the
  fractions of its own group (same material, different slice) but not across the two mixtures
  (different material, different TMT tube).

  @note @c MSstats_Mixture is an ordinary factor -- its name contains no "replicate" -- so it
        counts towards the condition. This design therefore has four OpenMS conditions
        (control/treated x mixture 1/2), not two, and the mixtures never share one. That affects
        only OpenMS' own condition grouping; the column is forwarded to MSstats unchanged.

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

  @note MS1-labeled designs have no MSstats route in OpenMS: @ref TOPP_MSstatsConverter
        @p -method @p LFQ refuses a design with more than one label and @p -method @p ISO requires
        @c MSstats_Mixture, and @ref TOPP_ProteomicsLFQ refuses multi-label designs outright.
        @ref TOPP_ProteinQuantifier consumes them. The MSstats columns above are shown only for
        consistency with the other examples.

  @subsection ExperimentalDesign_Ex7 7. The same design in two-table format

  Example 3 with the metadata factored out. The blank line separating the sections is mandatory.

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

  Every column of the sample section is a <b>factor</b>. Names are free-form; OpenMS only compares
  values, so @c Genotype, @c Timepoint or @c Dose behave like the MSstats columns below.

  A <b>condition</b> is a unique combination of the values of all factors except @c Sample and
  except every factor whose column name contains @c "replicate" or @c "Replicate". Samples that
  agree on all remaining factors form one condition. For a design whose sample section has
  factor columns, this is what getConditionToSampleMapping(), getSampleToConditionMapping(),
  getPathLabelToConditionMapping() and getConditionToPathLabelVector() return. For a factor-less
  section -- as built by fromConsensusMap() and fromIdentifications() -- they currently disagree:
  getSampleToConditionMapping() gives every sample its own condition, while the other three
  collapse all samples into one. @ref TOPP_Epifany uses getConditionToPathLabelVector() to merge
  ID runs per condition when a design is supplied; @ref TOPP_ProteomicsLFQ does not -- it merges
  every ID run study-wide regardless of condition.

  Condition numbers are the lexicographic rank of the factor-value tuple, so condition 0 is
  whichever value sorts first, not a reference level.

  @warning The replicate rule matches on the column NAME only. @c MSstats_BioReplicate and
           @c Technical_Replicate are excluded from the condition automatically; @c Donor or
           @c Rep are not, and every distinct value in them creates its own condition. This
           catches @c MSstats_Mixture too (see example 5). Name replicate columns accordingly.

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

  Checked when the file is loaded; a design that breaks one of these is rejected:

  - Column headers are matched exactly and are case-sensitive.
  - Fraction groups must be integers, must start at 1 and must be consecutive -- no gaps, no
    renumbering per condition.
  - Each <tt>(Fraction_Group, Fraction, Label)</tt> triple may occur only once.
  - Each <tt>(Spectra_Filepath, Label)</tt> pair may occur only once.
  - In a design with a single distinct @c Label value (normally label-free), each
    <tt>(Fraction_Group, Label)</tt> pair must map to exactly one sample. Two samples in one such
    fraction group means the labels are missing or the grouping is wrong.
  - One-table format only: @c Sample becomes mandatory as soon as any @c Label is greater than 1
    -- OpenMS cannot guess which channel is which material. Without a @c Sample column, the
    sample name defaults to the @c Fraction_Group value, which makes every fraction group its own
    sample.
  - Two-table format only: every @c Sample value used in the file section must exist in the
    sample section, otherwise loading fails with a bare @c std::out_of_range. Omitting the
    @c Sample column from a two-table file section therefore fails as well, unless the sample
    section literally contains rows named <tt>"Fraction group 1"</tt>,
    <tt>"Fraction group 2"</tt>, ...
  - A relative @c Spectra_Filepath is resolved first against the directory of the design file,
    then against the current working directory; if neither exists the string is kept as written.
    Tools that need the spectra present (@c require_spectra_files) fail instead.
  - A misspelled column name is not an error in the one-table format -- it silently becomes a
    sample-metadata factor. A mistyped @c Sample is the dangerous case: sample names then fall
    back to the fraction-group value, so reused samples quietly become distinct ones and the
    misspelled column also splits conditions. The two-table file section rejects unknown headers.
  - A design without data rows still gets its column headers checked; only the row-level rules
    above are skipped. A completely empty file loads as an empty design with no checks at all.

  @section ExperimentalDesign_SDRF Relation to SDRF-Proteomics

  <a href="https://github.com/bigbio/proteomics-sample-metadata">SDRF-Proteomics</a> is the
  HUPO-PSI sample-and-data-relationship format used by PRIDE and quantms. It covers more than an
  OpenMS design (organism, disease, instrument, search settings) and can generate one:
  <tt>parse_sdrf convert-openms -s sdrf.tsv</tt> from
  <a href="https://github.com/bigbio/sdrf-pipelines">sdrf-pipelines</a> writes the file described
  here. The terms correspond as follows -- for the exact conversion, which varies between
  releases, consult sdrf-pipelines itself:

  | OpenMS design           | SDRF-Proteomics                                              | Notes                                                                                     |
  |-------------------------|--------------------------------------------------------------|-------------------------------------------------------------------------------------------|
  | @c Spectra_Filepath     | <tt>comment[data file]</tt>                                  | copied verbatim; the converter can optionally rewrite the extension, e.g. <tt>.raw</tt> to <tt>.mzML</tt> |
  | @c Fraction             | <tt>comment[fraction identifier]</tt>                        | @c 1 when absent or "not available"                                                        |
  | @c Label                | <tt>comment[label]</tt>                                      | the ontology term becomes its 1-based index within the plex inferred for that file, so the same reagent name can map to different numbers in different plexes |
  | @c Fraction_Group       | <tt>source name</tt> + <tt>comment[technical replicate]</tt> | assigned per MS <i>file</i>, then renumbered to stay consecutive. Label-free data gets one group per (source, technical replicate) -- which is why re-injections get their own group; in a labeled design all channels of a file necessarily share its group |
  | @c Sample               | <tt>source name</tt>                                         | a source name ending in "sample N" contributes N, otherwise names are numbered in order of appearance |
  | @c MSstats_Condition    | the <tt>factor value[...]</tt> columns, joined with a pipe    | falls back to the non-redundant <tt>characteristics[...]</tt> columns, then to the source name   |
  | @c MSstats_BioReplicate | derived from <tt>source name</tt>                            | the source's sample number, so equal to @c Sample in a uniform design; <i>not</i> copied from <tt>characteristics[biological replicate]</tt> |
  | @c MSstats_Mixture      | derived per labeled file and sample                          | isobaric designs only (TMT/iTRAQ); SILAC gets no mixture column                             |

  Two SDRF terms do not map one-to-one. <tt>comment[technical replicate]</tt> is folded into
  @c Fraction_Group, since OpenMS expresses technical replication as "same sample, new fraction
  group" rather than as a column. And SDRF's <tt>source name</tt> (the material) vs.
  <tt>assay name</tt> (the measurement) corresponds conceptually to @c Sample vs. the file/label
  pair, though the OpenMS converter reads only the former.

  @warning Always check a converted design before using it, in particular for labeled data. The
           conversion is lossy in places and has had bugs -- sdrf-pipelines 0.1.6, for instance,
           writes @c Label @c 1 for every SILAC channel, which OpenMS then rejects as a duplicate
           <tt>(Spectra_Filepath, Label)</tt> pair.

  @section ExperimentalDesign_QPX Relation to QPX (quantms.io)

  QPX ("Quantitative Proteomics eXchange") is the Parquet exchange format written via the
  @c out_qpx option of @ref TOPP_ProteomicsLFQ, @ref TOPP_IsobaricWorkflow and
  @ref TOPP_ProteinQuantifier. OpenMS writes three views -- @c quantms.feature.parquet,
  @c quantms.psm.parquet and @c quantms.pg.parquet. It does <b>not</b> write @c run.parquet or
  @c sample.parquet; those are generated from the SDRF by the quantms converter.

  What the design controls in the views OpenMS writes:

  | OpenMS design                                 | QPX                                | Notes                                                                                   |
  |-----------------------------------------------|------------------------------------|-----------------------------------------------------------------------------------------|
  | basename of @c Spectra_Filepath, no extension | @c run_file_name                   | primary-key component of the psm and feature views; the pg view instead carries the run names in @c grouped_runs. See ArrowIOHelpers::qpxRunFileName() |
  | @c Fraction_Group                             | @c grouped_runs of the pg view     | QPX has no fraction-group number: the group is represented by listing its run names      |
  | @c Label                                      | <tt>intensities[].label</tt> (feature) and the scalar @c label (pg) | as a canonical channel token (@c "LFQ", @c "TMT126", <tt>"SILAC heavy"</tt>), not as the number; see ArrowIOHelpers::qpxIntensityLabels() |

  Those labels and run names are join keys against @c run.parquet, so the design also has to
  agree with the views OpenMS does not write: @c Fraction with @c run.fraction, @c Sample with
  <tt>run.samples[].sample_accession</tt> (keyed by the SDRF source name), and
  @c MSstats_BioReplicate with <tt>run.samples[].biological_replicate</tt>. QPX's
  <tt>run.samples[].technical_replicate</tt> has no counterpart in an OpenMS design, because
  technical replication is expressed as a second fraction group.

  @note The pg exporter requires a fraction group to be <b>rectangular</b>: every label the group
        publishes must exist in every one of its runs, so that <tt>(any run of the group, label)</tt>
        resolves to one sample. Ragged designs, a label meaning two samples within one group, and
        designs that put one run into several fraction groups are refused.

  @section ExperimentalDesign_Pitfalls Common pitfalls

  - Renumbering fraction groups per condition. They are global (1, 2, 3, ... across the whole
    design); restarting at 1 for the second condition fails to load.
  - Confusing fractions with replicates. Fractions are aggregated (one quantity from several
    files); replicates are kept apart (several quantities). Compare examples 2 and 3.
  - Writing one row per labeled file instead of one row per channel. Quantification then aborts:
    PeptideAndProteinQuant refuses a consensus column whose <tt>(file, channel)</tt> has no design
    row.
  - Naming a replicate column without "replicate" in it. It then counts towards the condition and
    splits each group of replicates into separate conditions.
  - Reusing a biological replicate id across conditions unintentionally. That declares a paired
    design; @ref TOPP_MSstatsConverter warns about it.
  - Saving as comma-separated. The format is TAB-separated.

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
      /// @brief Fraction group id: the rows combined into one quantity
      ///
      /// 1-based and consecutive over the whole design. A second measurement of the same
      /// material is a new group carrying the same @c sample.
      unsigned fraction_group = 1;

      /// @brief Fraction 1..m, mandatory, 1 if not set
      ///
      /// 1 everywhere for unfractionated data. The same number in two fraction groups denotes
      /// corresponding fractions.
      unsigned fraction = 1;

      /// @brief File name, mandatory
      std::string path = "UNKNOWN_FILE";

      /// @brief 1-based channel position within @c path (1 for label-free)
      ///
      /// E.g. 1..10 for TMT10plex in canonical channel order 126, 127N, ..., 131. A position,
      /// not a reagent name.
      unsigned label = 1;

      /// @brief Zero-based ROW INDEX into the sample section; allows grouping by sample
      ///
      /// This is not the @c Sample column value -- that is @c sample_name.
      unsigned sample = 0;

      /// @brief The @c Sample column value, i.e. the sample's NAME
      ///
      /// An arbitrary string such as "1" or "BSA1", used to look the sample row up by name
      /// rather than by index.
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
      /// @brief Get set of all factors (column names) defined for the sample section
      ///
      /// For a design parsed from a file this includes the @c Sample column itself, which the
      /// mapping functions of ExperimentalDesign drop explicitly before comparing rows. A section
      /// built with addSample() has no columns at all and returns an empty set.
      std::set< std::string > getFactors() const;

      /// Checks whether the sample section has a row for a sample name
      bool hasSample(const std::string& sample) const;

      /// Checks whether Sample Section has a specific factor (i.e. column name)
      bool hasFactor(const std::string &factor) const;

      /// Returns value of factor for given sample NAME and factor name
      std::string getFactorValue(const std::string& sample_name, const std::string &factor) const;

      /// Returns value of factor for given sample ROW INDEX (zero-based) and factor name
      std::string getFactorValue(unsigned sample_idx, const std::string &factor) const;

      /// @brief Returns column index of factor
      ///
      /// This is a storage detail, not a property of the design: a one-table design orders the
      /// sample columns alphabetically, a two-table design keeps the order of the file. Address
      /// factors by name unless the layout itself matters.
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

      @note If no factors are provided in the experimental design, each sample forms 
      its own unique condition (N samples = N distinct conditions).
      
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

    /// @brief Number of samples measured (= number of rows in the sample section)
    ///
    /// NOT the highest sample index: MSFileSectionEntry::sample is zero-based, so the highest
    /// index is one less than this.
    /// @return the number of rows in the sample section
    unsigned getNumberOfSamples() const;

    /// @brief Number of distinct fraction indices used anywhere in the design
    ///
    /// Meaningful only when all fraction groups use the same fraction numbering, which this does
    /// not check. sameNrOfMSFilesPerFraction() is a necessary but not sufficient sanity check: it
    /// only verifies that every fraction index occurs in equally many rows.
    /// @return the number of distinct fraction indices
    unsigned getNumberOfFractions() const;

    /// @brief Highest label index used anywhere in the design
    ///
    /// The plex size for a well-formed design, 1 for label-free. Contiguous labels are not
    /// enforced.
    /// @return the maximum label index, or 0 if the design is empty
    unsigned getNumberOfLabels() const;

    /// @brief Number of distinct MS file paths in the design
    ///
    /// A multiplexed file contributes one file but several rows.
    /// @return the number of distinct paths
    unsigned getNumberOfMSFiles() const;

    /// @brief Number of distinct fraction groups
    ///
    /// Allows to group fraction ids and source files.
    /// @return the number of distinct fraction groups
    unsigned getNumberOfFractionGroups() const;

    /// @brief Sample quantified in a given fraction group and label
    ///
    /// @return the zero-based sample row index quantified in @p fraction_group at @p label
    /// @throws Exception::ElementNotFound if the design has no such combination
    unsigned getSample(unsigned fraction_group, unsigned label = 1);

    /// @brief Whether the design is fractionated
    ///
    /// @return true if more than one distinct fraction index is used in the design
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
      @c sample_name as meta values on every column header whose <tt>(basename, label)</tt> pair
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

