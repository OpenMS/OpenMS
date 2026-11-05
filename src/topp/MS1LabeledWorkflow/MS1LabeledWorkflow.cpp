// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "MS1LabeledFAIMS.h"
#include "MS1LabeledRatioQuantifier.h"
#include "MS1LabeledSpectra.h"

#include <OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/BayesianProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/ConsensusMapMergerAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>
#include <OpenMS/ANALYSIS/QUANTITATION/DDAWorkflowCommons.h>
#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>
#include <OpenMS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>
#include <OpenMS/FEATUREFINDER/MultiplexResolverAlgorithm.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/FORMAT/ConsensusMapArrowExport.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/ProteinGroupArrowExport.h>
#include <OpenMS/FORMAT/QPXCollectionExport.h>
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/FORMAT/QPXIdentity.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/METADATA/MS1LabelState.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/SYSTEM/File.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_MS1LabeledWorkflow MS1LabeledWorkflow

@brief Complete quantification workflow for MS1-labeled (SILAC, Dimethyl, ...) LC-MS/MS experiments.

<CENTER>
  <table>
    <tr>
      <th ALIGN = "center"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=3> &rarr; MS1LabeledWorkflow &rarr;</td>
      <th ALIGN = "center"> pot. successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PercolatorAdapter </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_MzTabExporter </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter </td>
    </tr>
  </table>
</CENTER>

This tool runs the complete quantification of an experiment whose samples were labeled before the LC-MS
measurement so that the light and heavy forms of a peptide appear as separate MS1 features with a fixed mass
shift: SILAC (Lys4/Lys6/Lys8, Arg6/Arg10), Dimethyl and ICPL labeling, in duplex, triplex or higher plex.
It is the MS1-labeling counterpart of @ref TOPP_ProteomicsLFQ (label-free) and @ref TOPP_IsobaricWorkflow
(TMT / iTRAQ reporter ions) and combines the standalone tools @ref TOPP_FeatureFinderMultiplex,
@ref TOPP_IDMapper, @ref TOPP_IDConflictResolver, @ref TOPP_MultiplexResolver, @ref TOPP_MapAlignerIdentification,
@ref TOPP_FeatureLinkerUnlabeledQT, @ref TOPP_ProteinInference and @ref TOPP_ProteinQuantifier into one run.

<b>Input</b>
  - Spectra in mzML or Thermo RAW format, one file per LC-MS run (@p in), with unique basenames. Native
    RAW input requires a build with WITH_THERMO_RAW and a .NET 8+ runtime. MS1 peaks are used for
    quantification; spectrum metadata are also read to recover identification FAIMS CVs. Profile and
    centroided data are both accepted (see @p algorithm:spectrum_type).
  - Identifications (@p ids), one file per spectra file in the same order, already filtered at PSM level
    (e.g. q-value < 0.01) and carrying Posterior Error Probability scores, e.g. produced with
    @ref TOPP_PercolatorAdapter (@p -score_type pep) or @ref TOPP_IDPosteriorErrorProbability, followed by
    @ref TOPP_IDFilter. Pass the results of one search engine rather than concatenating those of several.
  - The label specification (@p labels), in the syntax of @ref TOPP_FeatureFinderMultiplex, e.g.
    <tt>[][Lys8,Arg10]</tt> for SILAC, <tt>[][Lys4,Arg6][Lys8,Arg10]</tt> for triple SILAC,
    <tt>[Dimethyl0][Dimethyl6]</tt> for Dimethyl. Every bracket is one channel; the channels are numbered
    from 1 in this order and are the @c Label column of the experimental design. Channels must be in
    increasing mass order, and SILAC requires an explicit unlabelled first channel <tt>[]</tt>.
  - Optionally an experimental design (@p design) with the columns Fraction_Group, Fraction,
    Spectra_Filepath, Label and Sample (see @ref OpenMS::ExperimentalDesign "ExperimentalDesign"). One row per
    (file, channel). Without a design every file is an unfractionated fraction group and every
    (file, channel) is its own sample.
  - Optionally a protein database (@p fasta). The identifications are then re-indexed with
    @ref TOPP_PeptideIndexer, which annotates protein sequences (for coverage), the decoy status and the
    theoretical peptide uniqueness needed by <tt>-protein_quantification strictly_unique_peptides</tt>.

FAIMS multiplets, identification mapping, blacklist checks and linking are restricted to the same
compensation voltage. CVs share the physical run's RT alignment and sample columns; they contribute
separate multiplet evidence to peptide quantification. Identification CVs are recovered from the original
spectrum references (including the preceding acquisition CV when MS2 metadata omit it). In a multi-CV
run, an ID without a resolvable spectrum reference must carry a valid @c FAIMS_CV annotation; it cannot
be assigned using RT alone. A single-CV run allows an unambiguous fallback to that voltage. Missing
references are repaired only against MS2 spectra at the known CV, within 0.01 seconds of the ID's RT.

<b>Important:</b> the labels have to be part of the database search as (variable) modifications, e.g.
<tt>Label:13C(6)15N(2) (K)</tt> and <tt>Label:13C(6)15N(4) (R)</tt> for Lys8/Arg10. Otherwise the MS2 spectra of the
labeled channels stay unidentified and the observed mass shifts cannot be reconciled with the peptide
sequences, so most multiplets end up as conflicts. The tool checks the search parameters recorded in @p ids for
the modifications implied by @p labels and refuses to run if they are missing (use @p -force to proceed anyway).

<b>Workflow</b>
  -# Per run: detection of peptide multiplets in the MS1 data (@ref OpenMS::FeatureFinderMultiplexAlgorithm,
     parameters in @p algorithm and @p label_mass_shifts). By default (@p algorithm:knock_out) partial multiplets are
     detected as well: a peptide seen in the heavy channel only, or whose weaker channels do not pass the detection's
     isotope pattern and intensity filters. A partial multiplet that merely flanks a larger multiplet of the same
     peptide (its channels at the same m/z, eluting next to it) is merged into the larger one.
  -# Per run: mapping of the identifications onto the multiplets (@p id_mapping). An identification matches a
     multiplet when it lies within the elution region of one of its channels (the RT range of its mass traces, plus
     @p id_mapping:rt_tolerance) and its m/z matches the channel - by default the theoretical m/z of the identified
     peptidoform, which carries the label, since the precursor m/z can be that of another isotopic peak. An
     identification of a channel the detection did not find (too weak) joins a partial multiplet detected in the
     peptide's other channels, at the m/z the peptide has there; its own channel is then re-quantified. With
     @p join_channels (default), the channels of an identified peptide that the detection found as separate partial multiplets
     (no spectrum held its complete pattern) are joined into one multiplet, and re-measured against its most intense channel
     from paired peaks over their common elution region, because each partial multiplet was quantified over other spectra.
  -# Per run: reduction to one identification per multiplet - preferring identifications whose labels agree with
     the multiplet's mass shifts over better-scoring ones that do not - and consolidation of quantitative and
     sequence information (@ref OpenMS::MultiplexResolverAlgorithm, @p resolver): multiplets whose mass shifts
     contradict the labels found in the annotated sequence are removed from quantification (their
     identifications are kept for protein inference), partial multiplets are completed with dummy
     features (intensity 0 = absent, NaN = not quantifiable, i.e. another detected multiplet overlaps). With
     @p requantify (default), a missing channel is then measured at its expected position over the elution region
     of the multiplet, like MaxQuant's re-quantify, and reported as absent only if it shows no signal there (its first
     two isotopes are not seen together in at least two spectra), and otherwise as not quantifiable if its peaks there do
     not show its isotope pattern in at least two spectra. The detection records which feature reserved each peak: only
     peaks reserved by another multiplet (also one of the same peptide, which is quantified on its own) make such a channel
     not quantifiable without measuring it, and of the reserved peaks it uses only those of patterns the detection
     dropped, not the isotopes of the multiplet's other channels; a channel that only the peaks it may not use would give a
     value is not quantifiable either, and so is a weak one (below @p requantify_correlation_ratio times the reference) whose
     elution profile correlates with the reference's by less than @p requantify_min_correlation (another species at its
     position; at moderate ratios real re-quantified channels correlate less well and are not checked). Re-quantified values are less accurate at extreme ratios
     than detected ones, so the ratios count the evidences that rest on them (see <b>Ratios</b> below). With
     @p targeted_quantification, an identified peptide that has no multiplet in
     the run (the detection missed it, or no multiplet fits its identifications) is measured at the m/z of its channels over
     the elution peak at its identification: the identified channel is integrated (if its isotope pattern fits the peptide),
     the other channels are measured against it from paired peaks (a channel is absent only if no peak is seen at its m/z).
     Such a feature and its identifications carry the meta
     value @c quantification_method = @c targeted. Its values are less accurate at extreme ratios than those of a detected
     multiplet, so they only add ratios and never change one that is reported without them (see <b>Ratios</b> below),
     and one whose other channel is not quantifiable is left out of the abundances. They take no part in the alignment, and
     they are linked among themselves only, after the detected multiplets, so that alignment, linking and match between runs
     stay as they are without them. One exception, with @p match_between_runs: a targeted feature replaces the unidentified
     multiplet kept for match between runs that it is measured from, so that run gets the targeted value instead of the one
     that linking would have given the multiplet (and that value no longer enters the run's normalization factors).
     Multiplets without identification are dropped,
     unless @p match_between_runs is set (see below). Once the resolver has used the labels, every identification
     is reduced to the <em>peptide identity</em>: the label modifications of @p labels are removed from the sequence,
     because the label belongs to the channel, not to the peptide, and the light and the heavy spectrum of one
     peptide have to name one peptide for linking, match between runs, inference and quantification (the
     convention MaxQuant uses as well). The label state stays documented on every identification as the meta
     values @c MS1Label:labeled_sequence (the peptidoform as searched), @c MS1Label:removed_labels (e.g. @c Lys8, or @c none) and
     @c MS1Label:channel (the 1-based channel the spectrum belongs to, i.e. the @c Label of the experimental
     design). The values also sit on every quantified consensus feature, for the identification it is quantified
     under. mzTab reports them as <tt>opt_global_*</tt> columns of the peptide and PSM sections, the QPX feature
     and psm views as @c cv_params. PSM-level output describes the spectrum match and therefore reports the
     peptidoform as searched (mzTab PSM section, QPX psm view), feature-level output the peptide identity
     (see @ref OpenMS::MS1LabelState). The column headers describe every channel's labels in @c channel_description.
     A spectrum match that was mapped onto several multiplets stays on the one whose matched channel is closest to
     the identified peptidoform in m/z; distinct spectra of one peptide on distinct multiplets are all kept, and their channel values
     add up per peptide and charge in the quantification.
  -# Per fraction: retention time alignment of the runs (@p alignment, identification-based, aligned to the
     run that shares the most identifications with every other run) and linking of the multiplets across runs (@p linking); the channels of
     every run are kept as sub-features, so the linked map has one column per (run, channel). Fractions are
     linked separately and then combined column-wise, exactly like ProteomicsLFQ does; a fraction measured in a
     single run is passed through. With @p match_between_runs, unidentified complete multiplets take part in the
     linking and take over the identification of a multiplet at the same position in another run (the SILAC
     equivalent of ProteomicsLFQ's <tt>-targeted_only false</tt>); multiplets that stay unidentified are not
     quantified. Unidentified partial multiplets are not kept: only an identification tells which channels they are.
  -# Protein inference over all runs (@p protein_inference), protein (and optionally PSM/peptide) FDR
     filtering, and peptide and protein quantification (@p ProteinQuantification), where the fractions of a
     fraction group are aggregated according to the design.

<b>Ratios.</b> A labeled experiment measures its channels in one run, so its quantity is their <em>ratio</em>,
and the tool computes it the way MaxQuant does, as a median of ratios rather than as a ratio of aggregated
intensities (@p ratios, computed by @c MS1LabeledRatioQuantifier next to this tool):
  - per multiplet and run, the intensity of a channel over that of the reference channel (@p ratios:reference_channel,
    the light one by default). Only positive, finite channels take part: an absent (dummy, intensity 0) or
    not-quantifiable channel is no measurement of a ratio.
  - per peptide and fraction group, the median of its evidence ratios.
  - per protein group and fraction group, the median of its peptides' ratios, reported only from
    @p ratios:min_ratio_count peptides upwards (MaxQuant's "min. ratio count", 2 by default), with the number of
    contributing peptides next to it (and, with @p targeted_quantification, the number of those from targeted quantification).
    Every ratio is also reported divided by the median peptide ratio of its
    (fraction group, channel), i.e. normalized on the assumption that most peptides do not change.
  - evidences of @p targeted_quantification only stand in where nothing else was measured: a peptide ratio uses them only if
    the peptide has no other evidence of that channel in that fraction group, a protein group ratio only if its other peptide ratios do
    not reach @p ratios:min_ratio_count, and the normalization factors only if a (fraction group, channel) has no other
    peptide ratio.

The reference channel is reported with the ratio 1.0 it has by construction, wherever another channel was
measured against it, so that every annotation covers the complete set of channels.

The ratios are annotated on the consensus features (@c MS1Label:evidence_ratio*, @c MS1Label:peptide_ratio*) and on the protein
groups, so they reach the consensusXML and the mzTab peptide section (as <tt>opt_global_*</tt> columns). In the
QPX pg view, whose rows are one per (protein group, fraction group, channel), they are written as that row's
@c additional_intensities, named @c ratio and @c ratio_normalized under the row's own channel label; the number
of contributing peptides sits in @c cv_params as @c ratio_count, being a count rather than an intensity, with
@p targeted_quantification next to the number of those from targeted quantification (@c ratio_targeted_count,
@c MS1Label:peptide_ratio_targeted_count on the features, @c fraction_group_level_ratio_targeted_count on the protein
groups). Where @p requantify measured channels, the number of contributing ratios that rest on them is reported as well
(@c ratio_requantified_count, @c MS1Label:peptide_ratio_requantified_count: the evidences whose channel or reference channel
was re-quantified; @c fraction_group_level_ratio_requantified_count: the peptide ratios all of whose evidences were), so that
they can be filtered at extreme ratios. No
separate @ref TOPP_ProteinQuantifier run is needed for any of it.

Next to the ratios, the per-channel <em>abundances</em> are reported as before (mzTab peptide and protein
sections, QPX @c intensities): the summed peptide intensities per channel, like MaxQuant's Intensity columns.
Dividing two of those is a ratio of aggregates, a different statistic from the ratios above, which weights
peptides by their intensity. @p ProteinQuantification:consensus:normalize scales every assay to the overall
median, which for a labeled experiment forces the median channel ratio to 1; leave it off unless that is intended.

@p max_nr_labelled_aas is used for both the feature detection and the resolver: it is the maximum number of
labelled amino acids per peptide minus one, i.e. for tryptic SILAC the number of allowed missed cleavages.
It should agree with the missed-cleavage setting of the search.

<b>Output</b> (at least one required; each output is optional individually)
  - consensusXML with the linked multiplets, one column per (run, channel) (@p out_cxml)
  - mzTab with peptide and protein abundances per assay (@p out)
  - QPX Parquet collection (@p out_qpx): quantms.feature.parquet, quantms.psm.parquet, quantms.pg.parquet.
    The channels are reported with the canonical SDRF/QPX labels (<tt>SILAC light</tt>, <tt>SILAC medium</tt>,
    <tt>SILAC heavy</tt>, <tt>DIMETHYL0</tt>, ...). Labels outside this vocabulary (ICPL, Leu3, plain mass
    shifts) cannot be exported to QPX; the tool refuses @p out_qpx for them up front.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_MS1LabeledWorkflow.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_MS1LabeledWorkflow.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMS1LabeledWorkflow :
  public TOPPBase
{
public:
  TOPPMS1LabeledWorkflow() :
    TOPPBase("MS1LabeledWorkflow", "Quantification workflow for MS1-labeled (SILAC, Dimethyl, ...) LC-MS/MS experiments.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<file list>", StringList(),
      "Input: spectra files (mzML, or Thermo RAW with WITH_THERMO_RAW and .NET 8+), one per LC-MS run. "
      "Only MS1 peaks are used; profile and centroided data are accepted.");
    setValidFormats_("in", ListUtils::create<std::string>("mzML"
#ifdef WITH_THERMO_RAW
      ",raw"
#endif
    ));

    registerInputFileList_("ids", "<file list>", StringList(),
      "Identifications filtered at PSM level (e.g., q-value < 0.01), one per spectra file in the same order.\n"
      "The identifications must carry Posterior Error Probability scores (e.g. PercolatorAdapter with -score_type pep,\n"
      "or IDPosteriorErrorProbability) and the labels must have been searched as (variable) modifications.\n"
      "Pass the results of one search engine rather than concatenating those of several.");
    setValidFormats_("ids", {"idXML", "mzId", "idparquet"});

    registerInputFile_("design", "<file>", "",
      "Experimental design (Fraction_Group, Fraction, Spectra_Filepath, Label, Sample), one row per (file, channel).\n"
      "'Label' is the 1-based position of the channel in '-labels'. If not given, every file is an unfractionated\n"
      "fraction group and every (file, channel) is a separate sample.", false);
    setValidFormats_("design", {"tsv"});

    registerInputFile_("fasta", "<file>", "",
      "Protein database. If given, the identifications are re-indexed (PeptideIndexer): protein sequences (for coverage),\n"
      "decoy annotation and theoretical peptide uniqueness (needed by 'strictly_unique_peptides') are taken from it.", false);
    setValidFormats_("fasta", {"fasta"});

    registerStringOption_("labels", "<text>", "[][Lys8,Arg10]",
      "Labels used for labelling the samples, one bracket per channel. [...] specifies the labels for a single sample. For example\n\n"
      "[][Lys8,Arg10]        ... SILAC\n"
      "[][Lys4,Arg6][Lys8,Arg10]        ... triple-SILAC\n"
      "[Dimethyl0][Dimethyl6]        ... Dimethyl\n"
      "[Dimethyl0][Dimethyl4][Dimethyl8]        ... triple Dimethyl\n"
      "[ICPL0][ICPL4][ICPL6][ICPL10]        ... ICPL\n"
      "The channels are numbered from 1 in this order ('Label' column of the experimental design).", false);
    registerIntOption_("max_nr_labelled_aas", "<int>", 0,
      "Maximum number of labelled amino acids per peptide, minus one. Peptides with up to (this value + 1) labelled amino acids\n"
      "are considered by feature detection and resolver. For SILAC with trypsin digestion, this is the maximum number of missed cleavages.", false);
    setMinInt_("max_nr_labelled_aas", 0);

    registerOutputFile_("out", "<file>", "", "Optional output mzTab file. At least one output must be specified.", false, false);
    setValidFormats_("out", {"mzTab"});
    registerOutputFile_("out_cxml", "<file>", "", "Optional output consensusXML file. At least one output must be specified.", false, false);
    setValidFormats_("out_cxml", {"consensusXML"});
    registerOutputDir_("out_qpx", "<directory>", "",
      "Optional output directory for QPX Parquet files (quantms.feature.parquet, quantms.psm.parquet, quantms.pg.parquet). At least one output must be specified.", false, false);

    registerDoubleOption_("proteinFDR", "<threshold>", 0.05, "Protein FDR threshold (0.05=5%).", false);
    setMinFloat_("proteinFDR", 0.0);
    setMaxFloat_("proteinFDR", 1.0);
    registerStringOption_("picked_proteinFDR", "<choice>", "false", "Use a picked protein FDR?", false);
    setValidStrings_("picked_proteinFDR", {"true", "false"});
    registerDoubleOption_("psmFDR", "<threshold>", 1.0,
      "FDR threshold for sub-protein level (e.g. 0.05=5%). Use -FDR_type to choose the level. Cutoff is applied at the highest level.", false);
    setMinFloat_("psmFDR", 0.0);
    setMaxFloat_("psmFDR", 1.0);
    registerStringOption_("FDR_type", "<option>", "PSM", "Sub-protein FDR level. PSM, PSM+peptide (best PSM q-value).", false);
    setValidStrings_("FDR_type", {"PSM", "PSM+peptide"});

    registerStringOption_("protein_inference", "<option>", "aggregation",
      "Infer proteins:\n"
      "aggregation  = aggregates all peptide scores across a protein (using the best score) \n"
      "bayesian     = computes a posterior probability for every protein based on a Bayesian network.", false);
    setValidStrings_("protein_inference", {"aggregation", "bayesian"});
    registerStringOption_("protein_quantification", "<option>", "unique_peptides",
      "Quantify proteins based on:\n"
      "unique_peptides = use peptides mapping to single proteins or a group of indistinguishable proteins"
      "(according to the set of experimentally identified peptides).\n"
      "strictly_unique_peptides = use peptides mapping to a unique single protein only.\n"
      "shared_peptides = use shared peptides only for its best group (by inference score)", false, true);
    setValidStrings_("protein_quantification", {"unique_peptides", "strictly_unique_peptides", "shared_peptides"});

    registerStringOption_("match_between_runs", "<option>", "false",
      "true: keep multiplets without an identification, so that linking can hand them the identification of a multiplet at the\n"
      "same position in another run (the counterpart of ProteomicsLFQ's '-targeted_only false').\n"
      "false: only identified multiplets are quantified.\n"
      "Only complete unidentified multiplets are kept: the channels of a partial one (see 'algorithm:knock_out') are only known from its identification.", false);
    setValidStrings_("match_between_runs", {"true", "false"});

    registerStringOption_("requantify", "<option>", "true",
      "true: a channel of an identified multiplet that the detection did not find (see 'algorithm:knock_out') is re-quantified from the\n"
      "MS1 peaks at its expected position over the elution region of the multiplet, relative to the most intense detected channel (like\n"
      "MaxQuant's 're-quantify'). It is reported as absent (0) only if its first two isotopic traces are not seen together in at least\n"
      "two spectra there. Otherwise it is not quantifiable if those traces are seen only on peaks reserved for other features or channels,\n"
      "or if fewer than two spectra paired with the most intense detected channel show its isotope pattern (both traces, and no peak one\n"
      "isotope lower that exceeds 30 % of its monoisotopic peak, unless that position lies within twice 'algorithm:mz_tolerance' of isotope\n"
      "M+0 to M+n of another channel, with n + 1 the upper bound of 'algorithm:isotopes_per_peptide').\n"
      "The number of ratios that rest on re-quantified channels is reported ('MS1Label:peptide_ratio_requantified_count',\n"
      "'fraction_group_level_ratio_requantified_count', QPX 'ratio_requantified_count'); see also 'requantify_min_correlation'.\n"
      "false: such channels are reported as absent.", false);
    setValidStrings_("requantify", {"true", "false"});
    registerDoubleOption_("requantify_min_correlation", "<value>", 0.5,
      "Minimum Pearson correlation between the elution profile of a weak re-quantified channel (see 'requantify_correlation_ratio') and\n"
      "that of the reference channel over the multiplet's elution region (summed isotopic peaks per spectrum). Below it, the channel is\n"
      "not quantifiable: the channels of one peptide co-elute, so a signal that does not follow the reference is another species.\n"
      "-1: no such check.", false, true);
    setMinFloat_("requantify_min_correlation", -1.0);
    setMaxFloat_("requantify_min_correlation", 1.0);
    registerDoubleOption_("requantify_correlation_ratio", "<value>", 0.3,
      "'requantify_min_correlation' applies to the re-quantified channels whose intensity is below this fraction of the reference\n"
      "channel's: at moderate ratios the re-quantified channels are reliable, and real ones there correlate less well (the\n"
      "detection already missed them). A large value checks every re-quantified channel.", false, true);
    setMinFloat_("requantify_correlation_ratio", 0.0);

    registerStringOption_("join_channels", "<option>", "true",
      "true: when the detection found the channels of an identified peptide as separate partial multiplets (see 'algorithm:knock_out'),\n"
      "they are joined into one multiplet: partial multiplets at the m/z the peptide has in its missing channels (within 'algorithm:mz_tolerance'),\n"
      "of the same charge and CV, eluting at the identification and identified as the same peptide or not at all. The joined channels are\n"
      "re-measured against the most intense one from paired peaks over their common elution region, since the two partial multiplets were\n"
      "quantified over different spectra; a candidate with fewer than six spectra of paired peaks is not joined, nor one with a re-measured\n"
      "channel that is absent or not quantifiable (e.g. its peaks do not show its isotope pattern, see 'requantify'): its own channels, or the\n"
      "identified multiplet's channels if the candidate holds the most intense one. The remaining candidates (in triple SILAC, the medium and\n"
      "the heavy channel can be two) are then re-measured in the same way; the identified multiplet keeps its own values unless one of them is\n"
      "joined and holds the most intense channel.\n"
      "false: such partial multiplets are quantified separately (a missing channel that lies on another multiplet of the same peptide is\n"
      "then not quantifiable).", false);
    setValidStrings_("join_channels", {"true", "false"});

    registerStringOption_("targeted_quantification", "<option>", "true",
      "true: an identified peptide that has no multiplet in the run (none was detected, or its identifications fit none) is quantified from the MS1 peaks at\n"
      "the m/z it has in every channel: over the elution peak at the identification, the identified channel is integrated, and every other channel is\n"
      "measured against it from paired peaks (absent if no peak is seen at its m/z, not quantifiable if another signal occupies its position or its peaks do\n"
      "not show its isotope pattern in at least two spectra paired with the identified channel). The identified channel has to show its first two isotopes\n"
      "in at least three spectra and an isotope pattern that fits the peptide ('algorithm:averagine_similarity'). A lower channel's isotope in a channel's\n"
      "peak (unresolved: the same centroid, e.g. the light M+4 in a Dimethyl4 or Lys4 channel) is subtracted, estimated from that channel's monoisotopic\n"
      "peak; a channel whose peaks are mostly (20 % or more) such isotopes is not quantifiable. Such features and their identifications carry the meta value\n"
      "'quantification_method' = 'targeted'. Peptides quantified by a multiplet in the run, identifications on a channel of another peptide's multiplet,\n"
      "decoys and identifications whose labels fit no pattern are not measured; of identifications of different peptides on one signal, only the best one is.\n"
      "These values are less accurate at extreme ratios, so they only add ratios: a peptide or protein group ratio that the other evidences\n"
      "give (see 'ratios:min_ratio_count') does not change, and a feature whose other channel is not quantifiable is left out of the abundances.\n"
      "They take no part in the alignment and are linked among themselves only. The number of targeted ratios behind every peptide and protein\n"
      "group ratio is reported ('MS1Label:peptide_ratio_targeted_count', 'fraction_group_level_ratio_targeted_count', QPX 'ratio_targeted_count').\n"
      "Exception with 'match_between_runs': a targeted feature replaces the unidentified multiplet kept for match between runs that it is\n"
      "measured from, so that run gets the targeted value instead of the one linking would have given that multiplet.\n"
      "false: such identifications are kept for protein inference only.", false);
    setValidStrings_("targeted_quantification", {"true", "false"});

    // feature detection: the FeatureFinderMultiplex parameters, minus the two that are set at tool level
    Param ffm_defaults = FeatureFinderMultiplexAlgorithm().getDefaults();
    Param detection_defaults = ffm_defaults.copy("algorithm:", true);
    detection_defaults.remove("labels");
    detection_defaults.remove("max_nr_labelled_aas");
    detection_defaults.setValue("knock_out", "true",
      "Detect partial multiplets as well, e.g. a peptide seen in the heavy channel only, or whose light channel is too weak to be detected.\n"
      "Their missing channels are re-quantified (see 'requantify') or reported as absent. With 'false', only multiplets in which every\n"
      "channel was detected are quantified. Supported for up to four channels (more: always 'false').");
    Param mass_shift_defaults = ffm_defaults.copy("labels:", true);

    // resolver: tolerances only, labels and max_nr_labelled_aas are shared with the feature detection
    Param resolver_defaults = MultiplexResolverAlgorithm().getDefaults().copy("algorithm:", true);
    resolver_defaults.remove("labels");
    resolver_defaults.remove("max_nr_labelled_aas");
    resolver_defaults.remove("use_blacklist_owners"); // set from 'requantify'

    Param idm_defaults = IDMapper().getDefaults();
    // The identified peptidoform carries its labels, so its theoretical m/z is the monoisotopic m/z of its channel. The
    // precursor m/z is off by one or more isotopes whenever the instrument picked another peak of the envelope, which a
    // search with isotope error correction accepts.
    idm_defaults.setValue("mz_reference", "peptide",
      "Source of the m/z values of the identifications: 'peptide' = the theoretical m/z of the identified peptidoform (labels included),\n"
      "'precursor' = the measured precursor m/z (fails for spectra whose precursor is not the monoisotopic peak).");
    idm_defaults.setValue("rt_tolerance", idm_defaults.getValue("rt_tolerance"),
      "RT tolerance (in seconds) beyond the elution region of a multiplet: an identification matches a multiplet if it lies within the\n"
      "region of one of its channels' mass traces (plus the tolerances), not only near its apex.");
    idm_defaults.addTag("mz_reference", "advanced");
    idm_defaults.addTag("ignore_charge", "advanced");

    Param ma_defaults = MapAlignmentAlgorithmIdentification().getDefaults();
    ma_defaults.setValue("max_rt_shift", 0.1);
    ma_defaults.setValue("use_unassigned_peptides", "false");
    ma_defaults.setValue("use_feature_rt", "true");
    for (const auto& s : {"score_type", "score_cutoff", "min_score", "use_unassigned_peptides", "use_feature_rt", "use_adducts", "auto_reference", "auto_reference_min_points"})
    {
      ma_defaults.addTag(s, "advanced");
    }

    Param fl_defaults = FeatureGroupingAlgorithmQT().getDefaults();
    fl_defaults.setValue("distance_MZ:max_difference", 10.0);
    fl_defaults.setValue("distance_MZ:unit", "ppm");
    fl_defaults.setValue("distance_MZ:weight", 5.0);
    fl_defaults.setValue("distance_intensity:weight", 0.1);
    fl_defaults.setValue("use_identifications", "true");
    fl_defaults.remove("distance_RT:max_difference"); // estimated from the alignment
    for (const auto& s : {"distance_MZ:weight", "distance_intensity:weight", "use_identifications", "ignore_charge", "ignore_adduct"})
    {
      fl_defaults.addTag(s, "advanced");
    }

    Param pq_defaults = PeptideAndProteinQuant().getDefaults();
    // The reported quantity of a labeled experiment is the ratio (see the 'ratios' section), so the
    // abundances next to it are plain summed intensities, as MaxQuant's Intensity columns are, rather
    // than the median of the three most abundant peptides that a label-free run defaults to.
    pq_defaults.setValue("top:N", 0);
    pq_defaults.setValue("top:aggregate", "sum");
    pq_defaults.setValue("top:include_all", "true");
    pq_defaults.addTag("top:include_all", "advanced");

    Param ratio_defaults = MS1LabeledRatioQuantifier().getDefaults();

    Param combined;
    combined.insert("algorithm:", detection_defaults);
    combined.setSectionDescription("algorithm", "Parameters of the multiplet detection (FeatureFinderMultiplex)");
    combined.insert("label_mass_shifts:", mass_shift_defaults);
    combined.setSectionDescription("label_mass_shifts", "Mass shifts of all labels that can be used in '-labels'");
    combined.insert("resolver:", resolver_defaults);
    combined.setSectionDescription("resolver", "Parameters of the multiplet completion and quant/ID conflict resolution (MultiplexResolver)");
    combined.insert("id_mapping:", idm_defaults);
    combined.setSectionDescription("id_mapping", "Parameters for mapping the identifications onto the multiplets (IDMapper)");
    combined.insert("alignment:", ma_defaults);
    combined.setSectionDescription("alignment", "Parameters of the identification-based retention time alignment (MapAlignerIdentification)");
    combined.insert("linking:", fl_defaults);
    combined.setSectionDescription("linking", "Parameters for linking the multiplets across runs (FeatureLinkerUnlabeledQT)");
    combined.insert("ProteinQuantification:", pq_defaults);
    combined.setSectionDescription("ProteinQuantification", "Parameters of the peptide and protein abundances (ProteinQuantifier)");
    combined.insert("ratios:", ratio_defaults);
    combined.setSectionDescription("ratios", "Parameters of the channel ratios, the reported quantity of a labeled experiment");

    registerFullParam_(combined);
  }

  //-------------------------------------------------------------
  // identification input
  //-------------------------------------------------------------

  ExitCodes checkSingleRunPerID_(const vector<ProteinIdentification>& protein_ids, const std::string& id_file_abs_path) const
  {
    if (protein_ids.size() != 1)
    {
      OPENMS_LOG_FATAL_ERROR << "Exactly one protein identification run must be annotated in " << id_file_abs_path << endl;
      return ExitCodes::INCOMPATIBLE_INPUT_DATA;
    }

    StringList run_paths;
    protein_ids[0].getPrimaryMSRunPath(run_paths);
    if (run_paths.size() > 1)
    {
      OPENMS_LOG_FATAL_ERROR << "MS1LabeledWorkflow does not support merged ID runs. ID file: " << id_file_abs_path << endl;
      return ExitCodes::INCOMPATIBLE_INPUT_DATA;
    }
    if (run_paths.empty())
    {
      OPENMS_LOG_WARN << "Warning: No spectra file origin annotated in ID file. This can lead to errors or unexpected behaviour later: "
                      << id_file_abs_path << endl;
    }
    return EXECUTION_OK;
  }

  ExitCodes switchScoreType_(PeptideIdentificationList& peptide_ids, const std::string& id_file_abs_path) const
  {
    try
    {
      IDScoreSwitcherAlgorithm switcher;
      Size c = 0;
      switcher.switchToGeneralScoreType(peptide_ids, IDScoreSwitcherAlgorithm::ScoreType::PEP, c);
    }
    catch (Exception::MissingInformation&)
    {
      OPENMS_LOG_FATAL_ERROR << "MS1LabeledWorkflow expects a Posterior Error Probability score in all Peptide IDs. ID file: "
                             << id_file_abs_path << endl;
      return ExitCodes::INCOMPATIBLE_INPUT_DATA;
    }
    return EXECUTION_OK;
  }

  /**
    @brief The labels of '-labels' must have been searched as modifications, otherwise the labeled channels stay unidentified.

    Compares the modification names recorded in the search parameters of @p run with the PSI-MS names of the
    labels (e.g. "Label:13C(6)15N(2)" for Lys8). Fails unless '-force' is set; only warns when the search
    parameters record no modifications at all, because then nothing can be verified.
  */
  ExitCodes checkLabelModifications_(const ProteinIdentification& run, const std::string& id_file_abs_path) const
  {
    if (required_label_mods_.empty()) { return EXECUTION_OK; } // plain mass shifts: nothing to verify

    const ProteinIdentification::SearchParameters& sp = run.getSearchParameters();
    vector<std::string> mods = sp.fixed_modifications;
    mods.insert(mods.end(), sp.variable_modifications.begin(), sp.variable_modifications.end());
    if (mods.empty())
    {
      OPENMS_LOG_WARN << "Warning: " << id_file_abs_path << " records no search modifications, so it cannot be verified that the labels "
                      << "were searched as modifications. If they were not, the labeled channels are unidentified and most multiplets "
                      << "will be reported as quant/ID conflicts." << endl;
      return EXECUTION_OK;
    }

    // "Label:13C(6)15N(2) (K)" -> "Label:13C(6)15N(2)"; compared by name so that "Label:13C(6)" (Arg6) is
    // not accepted on the strength of "Label:13C(6)15N(2)" (Lys8)
    set<std::string> searched;
    for (const std::string& mod : mods)
    {
      const auto sep = mod.rfind(" (");
      searched.insert(sep == std::string::npos ? mod : mod.substr(0, sep));
    }

    vector<std::string> missing;
    for (const auto& [short_name, long_name] : required_label_mods_)
    {
      if (!searched.contains(long_name))
      {
        missing.push_back(short_name + " = " + long_name);
      }
    }
    if (missing.empty()) { return EXECUTION_OK; }

    const std::string message = "The search parameters of " + id_file_abs_path + " do not contain the modification(s) implied by '-labels "
      + labels_ + "': " + ListUtils::concatenate(missing, ", ") + ". Without them the MS2 spectra of the labeled channels stay "
      "unidentified and the observed mass shifts cannot be reconciled with the peptide sequences, so most multiplets end up as "
      "conflicts. Search with the labels as variable modifications (e.g. 'Label:13C(6)15N(2) (K)' and 'Label:13C(6)15N(4) (R)' "
      "for Lys8/Arg10).";
    if (getFlag_("force"))
    {
      OPENMS_LOG_WARN << "Warning: " << message << " Continuing because of '-force'." << endl;
      return EXECUTION_OK;
    }
    OPENMS_LOG_FATAL_ERROR << message << " Use '-force' to run anyway." << endl;
    return ExitCodes::INCOMPATIBLE_INPUT_DATA;
  }

  /**
    @brief Reconcile one input file's decoy affix into the single study-wide one (see ProteomicsLFQ).
  */
  void recordDecoyAffix_(const std::string& decoy_string, bool decoy_prefix, bool decoy_inferred, const std::string& source_file)
  {
    if (!decoy_inferred) { return; } // no FASTA given, nothing was inferred from this file

    if (!picked_decoy_seen_)
    {
      picked_decoy_string_ = decoy_string;
      picked_decoy_prefix_ = decoy_prefix;
      picked_decoy_seen_ = true;
    }
    else if (decoy_string != picked_decoy_string_ || decoy_prefix != picked_decoy_prefix_)
    {
      OPENMS_LOG_WARN << "Warning: decoy affix inferred from " << source_file << " ('" << decoy_string << "', "
                      << (decoy_prefix ? "prefix" : "suffix") << ") differs from the one inferred from earlier input ('"
                      << picked_decoy_string_ << "', " << (picked_decoy_prefix_ ? "prefix" : "suffix")
                      << "). Keeping the first. Picked protein FDR assumes one affix for the whole study - check that all inputs "
                         "were searched against the same database.\n";
    }
  }

  /// Reduce identifications this tool cannot tell apart to one per spectrum (see ProteomicsLFQ).
  void reduceToOnePerSpectrum_(PeptideIdentificationList& peptide_ids, const std::string& id_file_abs_path) const
  {
    const auto report = IDConflictResolverAlgorithm::reduceToOnePerSpectrum(peptide_ids);

    if (report.removed > 0)
    {
      OPENMS_LOG_WARN << "Warning: " << report.removed << " identification(s) in " << id_file_abs_path
                      << " repeat a (spectrum, peptidoform, charge) already claimed by another identification, e.g. "
                      << report.example << ".\nKept the best-scoring one of each group; the others would have counted the same "
                         "measurement more than once. To control this upstream, pass the results of one search engine rather "
                         "than concatenating those of several." << endl;
    }
    if (report.inconsistent_score_direction > 0)
    {
      OPENMS_LOG_WARN << "Warning: " << report.inconsistent_score_direction << " group(s) of repeated identifications in "
                      << id_file_abs_path << " disagree on whether a higher score is better, so no best one could be chosen. "
                         "They were left as they are and will be counted more than once." << endl;
    }
    if (report.without_spectrum_reference > 0)
    {
      OPENMS_LOG_INFO << report.without_spectrum_reference << " identification(s) in " << id_file_abs_path
                      << " carry no spectrum reference and were therefore not checked for repeats." << endl;
    }
    if (report.multiply_identified_spectra > 0)
    {
      OPENMS_LOG_INFO << report.multiply_identified_spectra << " spectra in " << id_file_abs_path
                      << " carry more than one identification naming different peptidoforms (a chimeric spectrum, or search "
                         "engines that disagree). Each is quantified separately." << endl;
    }
  }

  /// Load one run's identifications and bring them into the shape the rest of the workflow expects (see ProteomicsLFQ).
  /// Retain the MS1 peaks read alongside scan metadata, so native RAW input is decoded only once.
  ExitCodes loadAndCleanupIDFile_(
    const std::string& id_file_abs_path,
    const std::string& mz_file,
    const std::string& in_db,
    const Size fraction_group,
    const Size fraction,
    vector<ProteinIdentification>& protein_ids,
    PeptideIdentificationList& peptide_ids,
    set<std::string>& fixed_modifications,    // adds to
    set<std::string>& variable_modifications, // adds to
    MSExperiment& ms1)
  {
    const std::string mz_file_abs_path = File::absolutePath(mz_file);
    FileHandler().loadIdentifications(id_file_abs_path, protein_ids, peptide_ids,
                                      {FileTypes::IDXML, FileTypes::MZIDENTML, FileTypes::IDPARQUET}, log_type_);

    ExitCodes e = checkSingleRunPerID_(protein_ids, id_file_abs_path);
    if (e != EXECUTION_OK) { return e; }

    e = checkLabelModifications_(protein_ids[0], id_file_abs_path);
    if (e != EXECUTION_OK) { return e; }

    // re-index
    if (!in_db.empty())
    {
      PeptideIndexing indexer;
      Param param_pi = indexer.getParameters();
      param_pi.setValue("missing_decoy_action", "silent");
      param_pi.setValue("write_protein_sequence", "true");
      param_pi.setValue("write_protein_description", "true");
      indexer.setParameters(param_pi);

      FASTAContainer<TFI_File> fasta_db(in_db);
      const PeptideIndexing::ExitCodes indexer_exit = indexer.run(fasta_db, protein_ids, peptide_ids);

      recordDecoyAffix_(indexer.getDecoyString(), indexer.isPrefix(), true, id_file_abs_path);
      if ((indexer_exit != PeptideIndexing::ExitCodes::EXECUTION_OK) && (indexer_exit != PeptideIndexing::ExitCodes::PEPTIDE_IDS_EMPTY))
      {
        if (indexer_exit == PeptideIndexing::ExitCodes::DATABASE_EMPTY) { return INPUT_FILE_EMPTY; }
        else if (indexer_exit == PeptideIndexing::ExitCodes::UNEXPECTED_RESULT) { return UNEXPECTED_RESULT; }
        else { return UNKNOWN_ERROR; }
      }
    }

    e = switchScoreType_(peptide_ids, id_file_abs_path);
    if (e != EXECUTION_OK) { return e; }

    IDFilter::keepBestPeptideHits(peptide_ids, false); // strict = false

    // add to the (global) set of fixed and variable modifications
    const vector<std::string>& var_mods = protein_ids[0].getSearchParameters().variable_modifications;
    const vector<std::string>& fixed_mods = protein_ids[0].getSearchParameters().fixed_modifications;
    variable_modifications.insert(var_mods.begin(), var_mods.end());
    fixed_modifications.insert(fixed_mods.begin(), fixed_mods.end());

    // 'protein_references' carries the theoretical uniqueness that 'strictly_unique_peptides' filters on
    const bool keep_protein_references = getStringOption_("protein_quantification") == "strictly_unique_peptides";

    // drop meta values that are not needed downstream to free some space
    for (PeptideIdentification& pid : peptide_ids)
    {
      for (PeptideHit& ph : pid.getHits())
      {
        vector<std::string> keys;
        ph.getKeys(keys);
        for (const auto& k : keys)
        {
          if (!(StringUtils::hasSubstring(k, "_score")
                || StringUtils::hasSubstring(k, "q-value")
                || StringUtils::hasPrefix(k, "Luciphor_global_flr") // MzTab reports it as the MS:1002380 localization score
                || k == "target_decoy" // keep target_decoy information for QC
                || (keep_protein_references && k == "protein_references")))
          {
            ph.removeMetaValue(k);
          }
        }
      }
    }

    // check and reannotate the MS run the identifications refer to
    StringList id_msfile_ref;
    protein_ids[0].getPrimaryMSRunPath(id_msfile_ref);
    if (id_msfile_ref.empty())
    {
      OPENMS_LOG_WARN << "MS run path not set in ID file: " << id_file_abs_path << endl
                      << "Resetting reference to MS file provided at same input position." << endl;
    }
    else if (id_msfile_ref.size() == 1)
    {
      const std::string in_bn = File::stemName(mz_file_abs_path);
      const std::string id_primaryMSRun_bn = File::stemName(id_msfile_ref[0]);
      if (in_bn != id_primaryMSRun_bn)
      {
        OPENMS_LOG_WARN << "MS run path referenced from ID file does not match MS file at same input position: " << id_file_abs_path << endl
                        << "Resetting reference to MS file provided at same input position." << endl;
      }
    }
    else
    {
      OPENMS_LOG_WARN << "Multiple MS files referenced from ID file: " << id_file_abs_path << endl
                      << "Resetting reference to MS file provided at same input position." << endl;
    }
    protein_ids[0].setPrimaryMSRunPath(StringList{mz_file});
    protein_ids[0].setMetaValue("fraction_group", fraction_group);
    protein_ids[0].setMetaValue("fraction", fraction);

    // make the run identifier unique across input files (users split spectra and ID files before running the analysis)
    const std::string old_identifier = protein_ids[0].getIdentifier();
    const std::string new_identifier = old_identifier + "_" + StringUtils::toStr(fraction_group) + "F" + StringUtils::toStr(fraction);
    protein_ids[0].setIdentifier(new_identifier);
    for (PeptideIdentification& p : peptide_ids)
    {
      if (p.getIdentifier() == old_identifier) { p.setIdentifier(new_identifier); }
      else { OPENMS_LOG_WARN << "Peptide ID identifier found not present in the protein ID" << endl; }
    }

    // Resolve FAIMS CVs while the original spectrum references are still available. An RT-based
    // repair cannot establish a CV reliably in a run that cycles between several voltages.
    MSExperiment spectrum_metadata;
    MS1LabeledSpectra::load(mz_file_abs_path, spectrum_metadata, ms1, log_type_);
    MS1LabeledFAIMS::annotateCompensationVoltages(spectrum_metadata, peptide_ids);

    // reannotate spectrum references if missing
    bool missing_spec_ref = false;
    for (const PeptideIdentification& pid : peptide_ids)
    {
      if (!pid.metaValueExists(Constants::UserParam::SPECTRUM_REFERENCE)
          || pid.getMetaValue(Constants::UserParam::SPECTRUM_REFERENCE).toString().empty())
      {
        missing_spec_ref = true;
        break;
      }
    }
    if (missing_spec_ref)
    {
      OPENMS_LOG_WARN << "Warning: Identification file " << id_file_abs_path
                      << " contains IDs without meta value for the spectrum native id.\n"
                         "OpenMS will try to reannotate them by matching retention times between ID and spectra." << endl;
      MS1LabeledSpectra::addMissingSpectrumReferences(spectrum_metadata, peptide_ids);
    }

    // deliberately last, i.e. after the spectrum-reference repair above
    reduceToOnePerSpectrum_(peptide_ids, id_file_abs_path);

    return EXECUTION_OK;
  }

  //-------------------------------------------------------------
  // per run: multiplet detection, ID mapping, resolving
  //-------------------------------------------------------------

  /**
    @brief Detect the peptide multiplets of one run in its MS1 spectra (FeatureFinderMultiplex, incl. its FAIMS handling)

    @p features receives the channel features with their convex hulls (the handles of @p multiplets refer to them by unique
    id), @p blacklist the peaks the detected patterns consist of and, with re-quantification, joining of channels or targeted
    quantification, @p centroided the centroided MS1 peaks the detection worked on. Spectra and features carry the FAIMS CV
    they were acquired at.
  */
  ExitCodes detectMultiplets_(const std::string& mz_file, MSExperiment&& exp, ConsensusMap& multiplets, FeatureMap& features,
                              MSExperiment& blacklist, MSExperiment& centroided)
  {
    if (exp.empty())
    {
      OPENMS_LOG_FATAL_ERROR << "No MS1 spectra found in " << mz_file << endl;
      return INPUT_FILE_EMPTY;
    }

    for (const auto& spec : exp)
    {
      if (IMTypes::determineIMFormat(spec) == IMFormat::IM_PEAK)
      {
        OPENMS_LOG_ERROR << "Error: " << mz_file << " contains per-peak ion mobility data (IM_PEAK, "
                         << imPeakTypeToString(spec.getIMPeakType())
                         << ") which is not supported by the multiplet detection. Preprocess with IonMobilityBinning or PeakPickerIM first." << endl;
        return INCOMPATIBLE_INPUT_DATA;
      }
    }

    Param ffm_param;
    ffm_param.insert("algorithm:", getParam_().copy("algorithm:", true));
    ffm_param.setValue("algorithm:labels", labels_);
    ffm_param.setValue("algorithm:max_nr_labelled_aas", max_nr_labelled_aas_);
    ffm_param.setValue("algorithm:knock_out", knock_out_ ? "true" : "false");
    ffm_param.insert("labels:", getParam_().copy("label_mass_shifts:", true));
    writeDebug_("Parameters passed to FeatureFinderMultiplexAlgorithm", ffm_param, 3);

    // FAIMS data is processed per compensation voltage (a single NaN-keyed group for non-FAIMS data)
    auto faims_groups = IMDataConverter::splitByFAIMSCV(std::move(exp));
    const bool has_faims = faims_groups.size() > 1 || !std::isnan(faims_groups[0].first);
    if (has_faims)
    {
      OPENMS_LOG_INFO << "FAIMS data detected with " << faims_groups.size() << " compensation voltage(s)." << endl;
    }

    bool first_group = true;
    for (auto& [group_cv, faims_group] : faims_groups)
    {
      FeatureFinderMultiplexAlgorithm algorithm;
      algorithm.setParameters(ffm_param);
      algorithm.setLogType(log_type_);
      algorithm.run(faims_group, true);

      ConsensusMap& consensus_cv = algorithm.getConsensusMap();
      if (first_group)
      {
        multiplets.setColumnHeaders(consensus_cv.getColumnHeaders());
        multiplets.setExperimentType(consensus_cv.getExperimentType());
        first_group = false;
      }
      for (auto& cf : consensus_cv)
      {
        if (has_faims) { cf.setMetaValue(Constants::UserParam::FAIMS_CV, group_cv); }
        multiplets.push_back(cf);
      }
      for (auto& feature : algorithm.getFeatureMap())
      {
        if (has_faims) { feature.setMetaValue(Constants::UserParam::FAIMS_CV, group_cv); }
        features.push_back(std::move(feature));
      }
      for (auto& spec : algorithm.getBlacklist())
      {
        if (has_faims) { spec.setMetaValue(Constants::UserParam::FAIMS_CV, group_cv); }
        blacklist.addSpectrum(std::move(spec));
      }
      if (requantify_ || join_channels_ || targeted_quantification_)
      {
        for (auto& spec : algorithm.getCentroidedExperiment())
        {
          if (has_faims) { spec.setMetaValue(Constants::UserParam::FAIMS_CV, group_cv); }
          centroided.addSpectrum(std::move(spec));
        }
      }
    }
    multiplets.ensureUniqueId();
    multiplets.sortByPosition();
    // the resolver and the re-quantification look peaks up by RT range, which needs the spectra of all FAIMS groups in one order
    blacklist.sortSpectra();
    blacklist.updateRanges();
    centroided.sortSpectra();
    centroided.updateRanges();

    OPENMS_LOG_INFO << "Detected " << multiplets.size() << " peptide multiplet(s) in " << File::basename(mz_file) << "." << endl;
    return EXECUTION_OK;
  }

  /// The FAIMS CV a multiplet, feature or spectrum was acquired at (none for conventional data)
  static std::optional<double> compensationVoltage_(const MetaInfoInterface& value)
  {
    if (!value.metaValueExists(Constants::UserParam::FAIMS_CV)) { return std::nullopt; }
    return static_cast<double>(value.getMetaValue(Constants::UserParam::FAIMS_CV));
  }

  /// The absolute m/z tolerance of the multiplet detection at @p mz
  double detectionMZTolerance_(double mz) const
  {
    const double tolerance = getParam_().getValue("algorithm:mz_tolerance");
    return getParam_().getValue("algorithm:mz_unit").toString() == "ppm" ? mz * tolerance * 1e-6 : tolerance;
  }

  /**
    @brief Merge partial multiplets into the larger multiplet of the same peptide they flank, and describe where every channel elutes.

    With knock-outs, a peptide whose weakest channel is detected only around its apex shows up as a complete multiplet there and
    as partial multiplets (without that channel) in its flanks, where the remaining channels are still detected. Such a partial
    multiplet - every channel of it coincides in m/z with a channel of a larger multiplet of the same charge and CV, and the two
    elute within half of 'algorithm:rt_typical' of each other - is removed, and its elution region is added to that of the
    larger one. Its intensities are not added: the larger multiplet quantifies all its channels over the same scans, so its
    ratios are not biased by flanks that lack a channel. Identifications in the flanks map onto the larger multiplet instead.

    @param[in,out] multiplets The detected multiplets; merged flanks are removed
    @param[in] features The channel features of the detection
    @param[out] merged_into For the channels of every merged flank (unique id), the coinciding channel of the multiplet it was
                merged into (for a flank of a flank, of the outermost one)
    @return One feature per (multiplet, channel): the unique id, position and charge of the channel's handle, and one convex hull
            per mass trace of the channel, spanning the elution region of the multiplet (incl. merged flanks). Identifications
            are mapped onto these, and missing channels are re-quantified over the same region.
  */
  FeatureMap mergeFragments_(ConsensusMap& multiplets, const FeatureMap& features, std::unordered_map<UInt64, UInt64>& merged_into) const
  {
    merged_into.clear();
    std::unordered_map<UInt64, const Feature*> feature_of;
    for (const Feature& f : features) { feature_of.emplace(f.getUniqueId(), &f); }

    // elution region of every multiplet: the RT range of the mass traces of all its channels
    vector<std::pair<double, double>> region(multiplets.size());
    for (Size m = 0; m < multiplets.size(); ++m)
    {
      double lo = multiplets[m].getRT(), hi = multiplets[m].getRT();
      for (const FeatureHandle& handle : multiplets[m].getFeatures())
      {
        if (const auto it = feature_of.find(handle.getUniqueId()); it != feature_of.end() && !it->second->getConvexHulls().empty())
        {
          const DBoundingBox<2> box = it->second->getConvexHull().getBoundingBox();
          lo = std::min(lo, box.minX());
          hi = std::max(hi, box.maxX());
        }
      }
      region[m] = {lo, hi};
    }

    // which larger multiplet each partial one flanks
    const double gap = 0.5 * static_cast<double>(getParam_().getValue("algorithm:rt_typical"));
    vector<std::pair<double, Size>> handle_mzs; // (m/z of a handle, its multiplet), sorted
    for (Size m = 0; m < multiplets.size(); ++m)
    {
      for (const FeatureHandle& handle : multiplets[m].getFeatures()) { handle_mzs.emplace_back(handle.getMZ(), m); }
    }
    std::sort(handle_mzs.begin(), handle_mzs.end());
    auto coincides = [this](const FeatureHandle& a, const FeatureHandle& b)
    { return std::abs(a.getMZ() - b.getMZ()) <= detectionMZTolerance_(a.getMZ()); };

    const Size none = std::numeric_limits<Size>::max();
    vector<Size> target(multiplets.size(), none);
    for (Size a = 0; a < multiplets.size(); ++a)
    {
      const ConsensusFeature& partial = multiplets[a];
      if (partial.size() >= multiplicity_ || partial.getFeatures().empty()) { continue; }
      const FeatureHandle& first = *partial.getFeatures().begin();
      const double tolerance = detectionMZTolerance_(first.getMZ());
      auto it = std::lower_bound(handle_mzs.begin(), handle_mzs.end(), std::make_pair(first.getMZ() - tolerance, Size(0)));
      for (; it != handle_mzs.end() && it->first <= first.getMZ() + tolerance; ++it)
      {
        const Size b = it->second;
        const ConsensusFeature& larger = multiplets[b];
        if (larger.size() <= partial.size() || larger.getCharge() != partial.getCharge()
            || compensationVoltage_(larger) != compensationVoltage_(partial))
        {
          continue;
        }
        if (region[a].first > region[b].second + gap || region[a].second < region[b].first - gap) { continue; }
        const bool contained = std::all_of(partial.getFeatures().begin(), partial.getFeatures().end(), [&](const FeatureHandle& h)
          { return std::any_of(larger.getFeatures().begin(), larger.getFeatures().end(), [&](const FeatureHandle& g) { return coincides(h, g); }); });
        if (!contained) { continue; }
        // the most complete candidate, then the closest one
        auto distance = [&](Size m) { return std::abs(0.5 * (region[m].first + region[m].second) - 0.5 * (region[a].first + region[a].second)); };
        if (target[a] == none || larger.size() > multiplets[target[a]].size()
            || (larger.size() == multiplets[target[a]].size() && distance(b) < distance(target[a])))
        {
          target[a] = b;
        }
      }
    }

    // a flank of a flank belongs to the outermost multiplet (sizes strictly increase along the chain)
    Size merged = 0;
    for (Size a = 0; a < multiplets.size(); ++a)
    {
      if (target[a] == none) { continue; }
      Size b = target[a];
      while (target[b] != none) { b = target[b]; }
      region[b].first = std::min(region[b].first, region[a].first);
      region[b].second = std::max(region[b].second, region[a].second);
      ++merged;
    }

    // every channel of a flank becomes the channel of the larger multiplet it coincides with (the closest one)
    for (Size a = 0; a < multiplets.size(); ++a)
    {
      if (target[a] == none) { continue; }
      for (const FeatureHandle& h : multiplets[a].getFeatures())
      {
        const FeatureHandle* closest = nullptr;
        for (const FeatureHandle& g : multiplets[target[a]].getFeatures())
        {
          if (coincides(h, g) && (closest == nullptr || std::abs(g.getMZ() - h.getMZ()) < std::abs(closest->getMZ() - h.getMZ()))) { closest = &g; }
        }
        if (closest != nullptr && h.hasValidUniqueId() && closest->hasValidUniqueId()) { merged_into[h.getUniqueId()] = closest->getUniqueId(); }
      }
    }
    for (auto& [flank, into] : merged_into)
    {
      for (auto it = merged_into.find(into); it != merged_into.end(); it = merged_into.find(into)) { into = it->second; }
    }

    vector<ConsensusFeature> kept;
    kept.reserve(multiplets.size() - merged);
    FeatureMap channels;
    for (Size m = 0; m < multiplets.size(); ++m)
    {
      if (target[m] != none) { continue; }
      for (const FeatureHandle& handle : multiplets[m].getFeatures())
      {
        Feature channel;
        channel.setUniqueId(handle.getUniqueId());
        channel.setRT(handle.getRT());
        channel.setMZ(handle.getMZ());
        channel.setCharge(handle.getCharge());
        channel.setIntensity(handle.getIntensity());
        if (const auto cv = compensationVoltage_(multiplets[m])) { channel.setMetaValue(Constants::UserParam::FAIMS_CV, *cv); }
        // one hull per mass trace of the channel (the monoisotopic one first), spanning the elution region of the multiplet
        vector<std::pair<double, double>> traces;
        if (const auto it = feature_of.find(handle.getUniqueId()); it != feature_of.end())
        {
          for (const ConvexHull2D& hull : it->second->getConvexHulls())
          {
            const DBoundingBox<2> box = hull.getBoundingBox();
            traces.emplace_back(box.minY(), box.maxY());
          }
        }
        if (traces.empty()) { traces.emplace_back(handle.getMZ(), handle.getMZ()); }
        for (const auto& [mz_lo, mz_hi] : traces)
        {
          ConvexHull2D hull;
          hull.addPoint(DPosition<2>(region[m].first, mz_lo));
          hull.addPoint(DPosition<2>(region[m].first, mz_hi));
          hull.addPoint(DPosition<2>(region[m].second, mz_lo));
          hull.addPoint(DPosition<2>(region[m].second, mz_hi));
          channel.getConvexHulls().push_back(hull);
        }
        channels.push_back(std::move(channel));
      }
      kept.push_back(std::move(multiplets[m]));
    }
    multiplets.resize(0);
    for (auto& cf : kept) { multiplets.push_back(std::move(cf)); }
    for (auto& [index, header] : multiplets.getColumnHeaders()) { header.size = 0; }
    for (const auto& cf : multiplets)
    {
      for (const FeatureHandle& handle : cf.getFeatures()) { ++multiplets.getColumnHeaders()[handle.getMapIndex()].size; }
    }
    channels.updateRanges();

    if (merged > 0)
    {
      OPENMS_LOG_INFO << merged << " partial multiplet(s) flank a larger multiplet of the same peptide and were merged into it." << endl;
    }
    return channels;
  }

  /// Spectra with paired peaks that a joined channel needs (joinChannelsOfIdentifiedPeptides_). On real SILAC data (PXD022659,
  /// 487 joins) the re-measured ratio was within 1.5x of the mixing ratio for 60 % of the joins with 3-5 such spectra, for 92 %
  /// with 6-9 and for 98 % with 10 or more.
  static constexpr Size min_join_paired_spectra_ = 6;

  /**
    @brief Join the channels of an identified peptide that the detection found as separate partial multiplets.

    When no spectrum holds a peptide's complete pattern, the detection can report its channels as separate partial multiplets.
    mergeFragments_ does not merge them (neither is larger), and the resolver would complete each with a dummy on the other's
    peaks. So an identified partial multiplet takes up the partial multiplets at the m/z its peptide has in its missing channels
    (MultiplexResolverAlgorithm::expectedChannelMZs(), within the detection's m/z tolerance, which its own channels have to
    meet as well): of the same charge and CV, eluting at the identification (within 'id_mapping:rt_tolerance'), with all their
    channels among the missing ones and identifications (if any) of the same peptide only (the sequence without labels, at the
    same charge). The identifications have to be settled before (keepMultiAssignedOnClosest_): a spectrum match that the
    mapping attached to several multiplets would let each of them take up the partner.

    The partial multiplets were quantified over different spectra and numbers of isotopes, so their intensities are not
    comparable. The most intense channel keeps its intensity, and the channels of the other partial multiplets are re-measured
    against it from paired peaks over the union of the two channels' elution regions, with as many mass traces as both have
    (requantifyChannel_; of the reserved peaks, a channel uses only its own and those of no reported feature). A candidate is not
    joined if a re-measured channel has no value (see requantifyChannel_; e.g. its peaks do not show its isotope pattern) or fewer
    than min_join_paired_spectra_ spectra of paired peaks (counted whether they show the pattern or not). The re-measured channel
    can also be one of the identified multiplet's own, if the candidate holds the most intense channel; that candidate is then not
    joined (an identified one stays a multiplet of its own). The remaining candidates (in triple SILAC, the medium and the heavy
    channel can be two) are then re-measured against the most intense channel among them and the identified multiplet, so the
    identified multiplet's channels keep their own values only if no candidate is joined or they hold that channel. The peak one
    isotope below a re-measured channel is not taken for another envelope if it is an isotope of another channel of the identified
    multiplet or its candidates (lowerIsChannel_()); a channel that none of them has is not in that list. No elution correlation is
    required: on real data it did not separate good from bad joins, the m/z tolerance and the re-measurement did.

    A partner is joined to one partial multiplet only. Another identified partial multiplet of the same channels that the partner
    also spans (e.g. one elution that the detection split into two, in RT, or in m/z where its centroids scatter by more than the
    detection's m/z tolerance) keeps its own multiplet, in which the partner's channel is not quantifiable: the partner's peaks
    belong to a multiplet of the same peptide now. It adds no ratio then, but the summed channel abundances miss the partner
    channel's share of it (or, if the partner was the reference, count its share of the identified channels twice).
  */
  void joinChannelsOfIdentifiedPeptides_(ConsensusMap& multiplets, const FeatureMap& channels, const MultiplexResolverAlgorithm& resolver,
                                         const MSExperiment& centroided, const MSExperiment& blacklist) const
  {
    std::unordered_map<UInt64, Size> multiplet_of;
    for (Size m = 0; m < multiplets.size(); ++m)
    {
      for (const FeatureHandle& handle : multiplets[m].getFeatures()) { multiplet_of.emplace(handle.getUniqueId(), m); }
    }
    std::unordered_map<UInt64, const Feature*> channel_of;
    vector<std::pair<double, const Feature*>> by_mz;
    by_mz.reserve(channels.size());
    for (const Feature& channel : channels)
    {
      channel_of.emplace(channel.getUniqueId(), &channel);
      by_mz.emplace_back(channel.getMZ(), &channel);
    }
    std::sort(by_mz.begin(), by_mz.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

    const double rt_tolerance = getParam_().getValue("id_mapping:rt_tolerance");
    // the joined channels' own peaks are reserved by them; without owners they cannot be told from others', so none is skipped
    const bool owners = MultiplexResolverAlgorithm::hasBlacklistOwners(blacklist);
    auto detected = [](const ConsensusFeature& cf)
    { return static_cast<Size>(std::count_if(cf.getFeatures().begin(), cf.getFeatures().end(), [](const FeatureHandle& h) { return h.hasValidUniqueId(); })); };
    // the channel whose expected m/z @p mz lies at
    auto channel_at = [this](const vector<double>& mzs, double mz) -> std::optional<Size>
    {
      for (Size j = 0; j < mzs.size(); ++j)
      {
        if (std::abs(mz - mzs[j]) <= detectionMZTolerance_(mzs[j])) { return j; }
      }
      return std::nullopt;
    };
    // the expected channel m/z of an identification (empty if its labels fit no pattern)
    auto expected_mzs = [&](const PeptideIdentification& id, Size& own) -> vector<double>
    {
      if (id.getHits().empty() || id.getHits()[0].getCharge() == 0) { return {}; }
      return resolver.expectedChannelMZs(id.getHits()[0].getSequence(), id.getHits()[0].getCharge(), own);
    };

    const Size none = std::numeric_limits<Size>::max();
    vector<bool> taken(multiplets.size(), false);    // took part in a join
    vector<bool> absorbed(multiplets.size(), false); // joined to another multiplet
    std::set<std::pair<Size, Size>> rejected;        // pairs not joined: too few paired spectra, or a re-measured channel without a value
    Size n_joined = 0;
    for (Size m = 0; m < multiplets.size(); ++m)
    {
      ConsensusFeature& target = multiplets[m];
      if (taken[m] || target.getPeptideIdentifications().empty() || detected(target) >= multiplicity_) { continue; }

      // the best identification is the one the multiplet will be quantified under
      const PeptideIdentification* best = nullptr;
      for (const PeptideIdentification& id : target.getPeptideIdentifications())
      {
        if (id.getHits().empty() || id.getHits()[0].getCharge() == 0) { continue; }
        const double score = id.getHits()[0].getScore();
        if (best == nullptr || (id.isHigherScoreBetter() ? score > best->getHits()[0].getScore() : score < best->getHits()[0].getScore())) { best = &id; }
      }
      if (best == nullptr) { continue; }
      Size own = 0;
      const vector<double> mzs = expected_mzs(*best, own);
      if (mzs.empty()) { continue; }
      const int charge = best->getHits()[0].getCharge();
      const double rt = best->getRT();
      const AASequence peptide = unlabelled_(best->getHits()[0]);
      const std::optional<double> cv = compensationVoltage_(target);

      // every channel of the multiplet lies at the m/z of the peptide in one of its channels
      vector<bool> present(mzs.size(), false);
      bool fits = true;
      for (const FeatureHandle& handle : target.getFeatures())
      {
        const std::optional<Size> j = channel_at(mzs, handle.getMZ());
        if (!j) { fits = false; break; }
        present[*j] = true;
      }
      if (!fits) { continue; }

      // the same peptide (sequence and other modifications, in another channel) at the same charge
      auto same_peptide = [&](const PeptideIdentification& id)
      {
        Size id_own = 0;
        const vector<double> id_mzs = expected_mzs(id, id_own);
        if (id_mzs.size() != mzs.size() || id.getHits()[0].getCharge() != charge || !(unlabelled_(id.getHits()[0]) == peptide)) { return false; }
        for (Size j = 0; j < mzs.size(); ++j)
        {
          if (std::abs(id_mzs[j] - mzs[j]) > detectionMZTolerance_(mzs[j])) { return false; }
        }
        return true;
      };

      // per missing channel, the closest partial multiplet with a channel there
      vector<Size> partners;
      for (Size j = 0; j < mzs.size(); ++j)
      {
        if (present[j]) { continue; }
        Size partner = none;
        double partner_error = std::numeric_limits<double>::max();
        const double tolerance = detectionMZTolerance_(mzs[j]);
        auto it = std::lower_bound(by_mz.begin(), by_mz.end(), mzs[j] - tolerance, [](const auto& a, double mz) { return a.first < mz; });
        for (; it != by_mz.end() && it->first <= mzs[j] + tolerance; ++it)
        {
          const Feature& channel = *it->second;
          if (channel.getCharge() != charge || compensationVoltage_(channel) != cv) { continue; }
          const DBoundingBox<2> box = channel.getConvexHull().getBoundingBox();
          if (rt < box.minX() - rt_tolerance || rt > box.maxX() + rt_tolerance) { continue; }
          const auto found = multiplet_of.find(channel.getUniqueId());
          if (found == multiplet_of.end()) { continue; }
          const Size c = found->second;
          if (c == m || taken[c] || rejected.contains(std::minmax(m, c)) || std::find(partners.begin(), partners.end(), c) != partners.end()) { continue; }
          const ConsensusFeature& candidate = multiplets[c];
          if (detected(candidate) >= multiplicity_) { continue; }
          if (!std::all_of(candidate.getFeatures().begin(), candidate.getFeatures().end(), [&](const FeatureHandle& h)
                           { const std::optional<Size> k = channel_at(mzs, h.getMZ()); return k && !present[*k]; }))
          {
            continue;
          }
          if (!std::all_of(candidate.getPeptideIdentifications().begin(), candidate.getPeptideIdentifications().end(), same_peptide)) { continue; }
          const double error = std::abs(channel.getMZ() - mzs[j]);
          if (error < partner_error)
          {
            partner = c;
            partner_error = error;
          }
        }
        if (partner == none) { continue; }
        for (const FeatureHandle& handle : multiplets[partner].getFeatures()) { present[*channel_at(mzs, handle.getMZ())] = true; }
        partners.push_back(partner);
      }

      // The channels of the partial multiplets other than that of the most intense channel are re-measured against it; a
      // partner that cannot be (or against whose channel the multiplet's own cannot be) is not joined.
      std::unordered_map<UInt64, double> remeasured;
      while (!partners.empty())
      {
        vector<std::pair<FeatureHandle, Size>> handles; // with their partial multiplet: 0 = target, i + 1 = partners[i]
        for (Size g = 0; g <= partners.size(); ++g)
        {
          for (const FeatureHandle& handle : multiplets[g == 0 ? m : partners[g - 1]].getFeatures()) { handles.emplace_back(handle, g); }
        }
        vector<double> channel_mzs; // all channels of the multiplet and its partners
        for (const auto& [handle, group] : handles) { channel_mzs.push_back(handle.getMZ()); }
        UInt64 measured = UniqueIdInterface::INVALID; // the channel being re-measured
        auto skip_reserved = [&](const MSSpectrum& black, Size index)
        { return owners && !mayUseReserved_(MultiplexResolverAlgorithm::blacklistOwner(black, index), measured); };
        const auto reference = std::max_element(handles.begin(), handles.end(), [](const auto& a, const auto& b)
                                                { return a.first.getIntensity() < b.first.getIntensity(); });
        const Feature* reference_channel = channel_of.at(reference->first.getUniqueId());
        const DBoundingBox<2> reference_box = reference_channel->getConvexHull().getBoundingBox();

        Size failed = none;
        remeasured.clear();
        for (const auto& [handle, group] : handles)
        {
          if (group == reference->second) { continue; }
          measured = handle.getUniqueId();
          const Feature* channel = channel_of.at(handle.getUniqueId());
          const DBoundingBox<2> box = channel->getConvexHull().getBoundingBox();
          const Size traces = std::min(channel->getConvexHulls().size(), reference_channel->getConvexHulls().size());
          const Requantified_ result = requantifyChannel_(centroided, blacklist, skip_reserved, cv, std::min(box.minX(), reference_box.minX()),
                                                          std::max(box.maxX(), reference_box.maxX()), handle.getMZ(), handle.getCharge(), channel_mzs,
                                                          traces, reference->first.getMZ(), reference->first.getIntensity());
          if (!(result.intensity > 0.0) || result.paired < min_join_paired_spectra_)
          {
            failed = group;
            break;
          }
          remeasured[handle.getUniqueId()] = result.intensity;
        }
        if (failed == none) { break; }
        const Size drop = failed != 0 ? failed : reference->second; // a partner's group, never the target's
        rejected.insert(std::minmax(m, partners[drop - 1]));
        partners.erase(partners.begin() + static_cast<std::ptrdiff_t>(drop - 1));
      }
      if (partners.empty()) { continue; }

      // the partners' channels and identifications join the multiplet
      vector<FeatureHandle> handles(target.getFeatures().begin(), target.getFeatures().end());
      for (Size partner : partners)
      {
        handles.insert(handles.end(), multiplets[partner].getFeatures().begin(), multiplets[partner].getFeatures().end());
        for (PeptideIdentification& id : multiplets[partner].getPeptideIdentifications())
        {
          id.setMetaValue("map_index", 0); // the channel is settled by keepMultiAssignedOnClosest_
          target.getPeptideIdentifications().push_back(std::move(id));
        }
        multiplets[partner].getPeptideIdentifications().clear();
        taken[partner] = absorbed[partner] = true;
        ++n_joined;
      }
      taken[m] = true;
      for (FeatureHandle& handle : handles)
      {
        if (const auto it = remeasured.find(handle.getUniqueId()); it != remeasured.end()) { handle.setIntensity(it->second); }
      }
      // channels in m/z order, numbered from 0 like a detected pattern; the lightest is the anchor
      std::sort(handles.begin(), handles.end(), [](const FeatureHandle& a, const FeatureHandle& b) { return a.getMZ() < b.getMZ(); });
      ConsensusFeature::HandleSetType renumbered;
      for (Size i = 0; i < handles.size(); ++i)
      {
        handles[i].setMapIndex(i);
        renumbered.insert(handles[i]);
      }
      target.setFeatures(std::move(renumbered));
      target.setMZ(handles.front().getMZ());
      target.setIntensity(handles.front().getIntensity());
    }

    if (n_joined > 0)
    {
      vector<ConsensusFeature> kept;
      kept.reserve(multiplets.size() - n_joined);
      for (Size m = 0; m < multiplets.size(); ++m)
      {
        if (!absorbed[m]) { kept.push_back(std::move(multiplets[m])); }
      }
      multiplets.resize(0);
      for (auto& cf : kept) { multiplets.push_back(std::move(cf)); }
      for (auto& [index, header] : multiplets.getColumnHeaders()) { header.size = 0; }
      for (const auto& cf : multiplets)
      {
        for (const FeatureHandle& handle : cf.getFeatures()) { ++multiplets.getColumnHeaders()[handle.getMapIndex()].size; }
      }
    }
    if (n_joined + rejected.size() > 0)
    {
      OPENMS_LOG_INFO << n_joined << " partial multiplet(s) were found to be other channels of an identified partial multiplet and were "
                      << "joined to it, re-measured over common spectra; " << rejected.size() << " such candidate(s) were not joined (fewer than "
                      << min_join_paired_spectra_ << " spectra with paired peaks, or a re-measured channel without a value)." << endl;
    }
  }

  /**
    @brief Attach identifications of channels the detection did not find to the multiplet of their peptide.

    An identification maps onto a multiplet through its own channel. When that channel was too weak to be detected (e.g. the
    minor channel of a mixture, which DDA fragments all the same), the multiplet can still be there, detected in the other
    channels of the peptide. So an identification that stayed unassigned is looked up at the m/z its peptide has in the other
    channels (MultiplexResolverAlgorithm::expectedChannelMZs()): if a partial multiplet of the same charge and CV has a channel
    there that elutes at the identification (within the 'id_mapping' tolerances), the identification is attached to it, and its
    own channel is added to the multiplet as a feature of intensity 0 for the re-quantification to measure. The resolver then
    checks the result like any other identified multiplet.
  */
  void attachToUndetectedChannels_(ConsensusMap& multiplets, const FeatureMap& channels, const MultiplexResolverAlgorithm& resolver) const
  {
    std::unordered_map<UInt64, Size> multiplet_of;
    for (Size m = 0; m < multiplets.size(); ++m)
    {
      for (const FeatureHandle& handle : multiplets[m].getFeatures()) { multiplet_of.emplace(handle.getUniqueId(), m); }
    }
    vector<std::pair<double, const Feature*>> by_mz;
    by_mz.reserve(channels.size());
    for (const Feature& channel : channels) { by_mz.emplace_back(channel.getMZ(), &channel); }
    std::sort(by_mz.begin(), by_mz.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

    const double rt_tolerance = getParam_().getValue("id_mapping:rt_tolerance");
    const double mz_tolerance = getParam_().getValue("id_mapping:mz_tolerance");
    const bool ppm = getParam_().getValue("id_mapping:mz_measure").toString() == "ppm";
    auto tolerance_at = [&](double mz) { return ppm ? mz * mz_tolerance * 1e-6 : mz_tolerance; };

    PeptideIdentificationList unassigned;
    Size attached = 0;
    for (PeptideIdentification& id : multiplets.getUnassignedPeptideIdentifications())
    {
      if (id.getHits().empty() || id.getHits()[0].getCharge() == 0)
      {
        unassigned.push_back(std::move(id));
        continue;
      }
      const PeptideHit& hit = id.getHits()[0];
      Size own = 0;
      const vector<double> mzs = resolver.expectedChannelMZs(hit.getSequence(), hit.getCharge(), own);
      const std::optional<double> cv = compensationVoltage_(id);

      // the closest channel feature at the m/z of the peptide in another channel
      const Feature* best = nullptr;
      double best_error = std::numeric_limits<double>::max();
      for (Size j = 0; j < mzs.size(); ++j)
      {
        if (j == own) { continue; }
        const double tolerance = tolerance_at(mzs[j]);
        auto it = std::lower_bound(by_mz.begin(), by_mz.end(), mzs[j] - tolerance, [](const auto& a, double mz) { return a.first < mz; });
        for (; it != by_mz.end() && it->first <= mzs[j] + tolerance; ++it)
        {
          const Feature& channel = *it->second;
          if (channel.getCharge() != hit.getCharge() || compensationVoltage_(channel) != cv) { continue; }
          const DBoundingBox<2> box = channel.getConvexHull().getBoundingBox();
          if (id.getRT() < box.minX() - rt_tolerance || id.getRT() > box.maxX() + rt_tolerance) { continue; }
          // in a complete multiplet the identification's channel was detected and did not match it (placeholders added
          // below for earlier identifications have no unique id and do not count)
          const auto& handles = multiplets[multiplet_of.at(channel.getUniqueId())].getFeatures();
          if (static_cast<Size>(std::count_if(handles.begin(), handles.end(), [](const FeatureHandle& h) { return h.hasValidUniqueId(); })) >= multiplicity_)
          {
            continue;
          }
          const double error = std::abs(channel.getMZ() - mzs[j]);
          if (error < best_error)
          {
            best = &channel;
            best_error = error;
          }
        }
      }
      if (best == nullptr)
      {
        unassigned.push_back(std::move(id));
        continue;
      }

      // the identification's own channel joins the multiplet (unless an earlier identification added it already)
      ConsensusFeature& multiplet = multiplets[multiplet_of.at(best->getUniqueId())];
      const double own_mz = mzs[own];
      vector<FeatureHandle> handles(multiplet.getFeatures().begin(), multiplet.getFeatures().end());
      if (std::none_of(handles.begin(), handles.end(), [&](const FeatureHandle& h) { return std::abs(h.getMZ() - own_mz) <= tolerance_at(own_mz); }))
      {
        FeatureHandle placeholder;
        placeholder.setMZ(own_mz);
        placeholder.setRT(multiplet.getRT());
        placeholder.setCharge(hit.getCharge());
        placeholder.setIntensity(0.0);
        handles.push_back(placeholder);
        // channels in m/z order, numbered from 0 like a detected (partial) pattern; the lightest is the anchor
        std::sort(handles.begin(), handles.end(), [](const FeatureHandle& a, const FeatureHandle& b) { return a.getMZ() < b.getMZ(); });
        ConsensusFeature::HandleSetType renumbered;
        for (Size i = 0; i < handles.size(); ++i)
        {
          handles[i].setMapIndex(i);
          renumbered.insert(handles[i]);
        }
        multiplet.setFeatures(std::move(renumbered));
        multiplet.setMZ(handles.front().getMZ());
      }
      // the channel is settled by keepMultiAssignedOnClosest_
      id.setMetaValue("map_index", 0);
      multiplet.getPeptideIdentifications().push_back(std::move(id));
      ++attached;
    }
    multiplets.getUnassignedPeptideIdentifications() = std::move(unassigned);

    if (attached > 0)
    {
      OPENMS_LOG_INFO << attached << " identification(s) of channels the detection did not find were attached to the multiplet of their "
                      << "peptide, detected in another channel." << endl;
    }
  }

  /**
    @brief Attach every spectrum match to the closest channel, and to the closest multiplet only.

    The mapping attaches an identification to every multiplet within its tolerances and records the
    channel it matched as its 'map_index'. The resolver reads that index as the channel the spectrum
    belongs to, so it is set to the channel closest to the identification's m/z here: the theoretical
    m/z of the identified peptidoform, whose labels name its channel (the precursor m/z can be off by
    isotopes and then lie closer to another channel, e.g. Lys4 at charge 3). One spectrum is evidence
    for one multiplet, too (and a QPX PSM may reference one feature): it stays with the multiplet whose
    closest channel is closest to the identification in m/z, then in RT. This runs before the best-scoring
    match is chosen per multiplet, so that a multiplet does not first take a neighbour's
    better-scoring spectrum and then lose it, ending up without its own. Distinct spectra of the same
    peptide on distinct multiplets are not touched: both multiplets are quantified and their channel
    values add up per peptide.
  */
  void keepMultiAssignedOnClosest_(ConsensusMap& multiplets) const
  {
    struct Candidate
    {
      Size feature;
      Size id_index;
      double mz_ppm;
      double rt_diff;
    };
    std::map<std::string, vector<Candidate>> by_spectrum_match; // (spectrum, peptidoform, charge) -> where it sits

    for (Size f = 0; f < multiplets.size(); ++f)
    {
      ConsensusFeature& cf = multiplets[f];
      auto& ids = cf.getPeptideIdentifications();
      for (Size i = 0; i < ids.size(); ++i)
      {
        PeptideIdentification& id = ids[i];
        if (id.getHits().empty()) { continue; }

        // the channel closest to the identified peptidoform is the one the spectrum belongs to
        const PeptideHit& hit = id.getHits()[0];
        const double id_mz = hit.getCharge() != 0 ? hit.getSequence().getMZ(hit.getCharge()) : id.getMZ();
        double best_ppm = std::numeric_limits<double>::infinity(), best_rt = 0.0;
        Size best_map_index = 0;
        for (const FeatureHandle& handle : cf.getFeatures())
        {
          const double ppm = std::abs(handle.getMZ() - id_mz) / handle.getMZ() * 1e6;
          if (ppm < best_ppm)
          {
            best_ppm = ppm;
            best_rt = std::abs(handle.getRT() - id.getRT());
            best_map_index = handle.getMapIndex();
          }
        }
        if (!cf.getFeatures().empty()) { id.setMetaValue("map_index", best_map_index); }

        const std::string reference = id.getSpectrumReference();
        if (reference.empty()) { continue; }
        const std::string key = reference + '\t' + id.getHits()[0].getSequence().toString() + '\t' + StringUtils::toStr(id.getHits()[0].getCharge());
        by_spectrum_match[key].push_back({f, i, best_ppm, best_rt});
      }
    }

    vector<std::set<Size>> to_remove(multiplets.size()); // per multiplet: indices of the copies to drop
    Size shared = 0;
    for (auto& [key, candidates] : by_spectrum_match)
    {
      if (candidates.size() < 2) { continue; }
      ++shared;
      const auto closest = std::min_element(candidates.begin(), candidates.end(),
                                            [](const Candidate& a, const Candidate& b)
                                            { return std::tie(a.mz_ppm, a.rt_diff) < std::tie(b.mz_ppm, b.rt_diff); });
      for (auto it = candidates.begin(); it != candidates.end(); ++it)
      {
        if (it != closest) { to_remove[it->feature].insert(it->id_index); }
      }
    }
    for (Size f = 0; f < multiplets.size(); ++f)
    {
      if (to_remove[f].empty()) { continue; }
      auto& ids = multiplets[f].getPeptideIdentifications().getData();
      // highest index first, so that the remaining indices stay valid
      for (auto it = to_remove[f].rbegin(); it != to_remove[f].rend(); ++it)
      {
        ids.erase(ids.begin() + static_cast<std::ptrdiff_t>(*it));
      }
    }
    if (shared > 0)
    {
      OPENMS_LOG_INFO << shared << " spectrum match(es) were mapped onto several multiplets; each stays with the closest one." << endl;
    }
  }

  /**
    @brief Of several identifications on one multiplet, keep those whose labels agree with its mass shifts.

    The best-scoring identification is the one a multiplet is quantified under, and if its labels contradict the detected mass
    shifts, the multiplet is a quant/ID conflict - even when another identification on it fits (e.g. the best one comes from
    a co-isolated peptide). So the identifications that fit are preferred: those that do not are moved to the unassigned ones
    (they remain PSMs for protein inference). A multiplet that none of its identifications fits stays a conflict.
  */
  static void keepConsistentIdentifications_(ConsensusMap& multiplets, const MultiplexResolverAlgorithm& resolver)
  {
    Size moved = 0;
    for (auto& cf : multiplets)
    {
      PeptideIdentificationList& ids = cf.getPeptideIdentifications();
      if (ids.size() < 2) { continue; }
      const vector<bool> consistent = resolver.checkIdentifications(cf);
      if (std::none_of(consistent.begin(), consistent.end(), [](bool c) { return c; })) { continue; }
      if (std::all_of(consistent.begin(), consistent.end(), [](bool c) { return c; })) { continue; }
      PeptideIdentificationList kept;
      for (Size i = 0; i < ids.size(); ++i)
      {
        if (consistent[i])
        {
          kept.push_back(std::move(ids[i]));
          continue;
        }
        ids[i].removeMetaValue("map_index");
        multiplets.getUnassignedPeptideIdentifications().push_back(std::move(ids[i]));
        ++moved;
      }
      ids = std::move(kept);
    }
    if (moved > 0)
    {
      OPENMS_LOG_INFO << moved << " identification(s) contradict the labels of the multiplet they were mapped onto, which another "
                      << "identification on it fits; they are not used for quantification." << endl;
    }
  }

  /// May a channel use a peak that the detection reserved for @p owner (blacklist with owners)? Only if the channel reserved it itself
  /// (@p channel: its unique id, INVALID if it was not detected) or no reported feature did: another channel's peak at its position is
  /// that channel's isotope, and another multiplet quantifies its peaks itself (MultiplexResolverAlgorithm::isForeign()).
  static bool mayUseReserved_(UInt64 owner, UInt64 channel)
  {
    return owner == UniqueIdInterface::INVALID || owner == channel;
  }

  /// The blacklist spectrum of @p spectrum (same RT and CV), if any
  static const MSSpectrum* blacklistSpectrum_(const MSExperiment& blacklist, const MSSpectrum& spectrum, const std::optional<double>& cv)
  {
    for (auto b = blacklist.RTBegin(spectrum.getRT()); b != blacklist.end() && b->getRT() == spectrum.getRT(); ++b)
    {
      if (compensationVoltage_(*b) == cv) { return &*b; }
    }
    return nullptr;
  }

  /// A peak one isotope lower than a channel that exceeds this fraction of the channel's monoisotopic peak belongs to another envelope,
  /// whose isotopes the channel's peaks would be (as MultiplexFiltering::checkForSignificantPeak_)
  static constexpr double lower_isotope_fraction_ = 0.3;
  /// Spectra in which a channel measured against a reference needs its isotope pattern (first two isotopes, and no peak one isotope
  /// lower that exceeds lower_isotope_fraction_ of its monoisotopic peak) paired with the reference
  static constexpr Size min_pattern_spectra_ = 2;

  /// The largest isotope (M+n) the detection assigns to a channel: the maximum of 'algorithm:isotopes_per_peptide' minus one
  Size maxIsotope_() const
  {
    const std::string isotopes_per_peptide = getParam_().getValue("algorithm:isotopes_per_peptide").toString();
    return static_cast<Size>(std::max(StringUtils::toInt32(StringUtils::prefix(isotopes_per_peptide, ':')),
                                      StringUtils::toInt32(StringUtils::suffix(isotopes_per_peptide, ':'))) - 1);
  }

  /**
    @brief Is @p lower_mz, one isotope below channel @p self of the peptide's channels @p mzs (charge @p charge), an isotope of another of them?

    True if it lies within twice the detection's m/z tolerance (the window in which MultiplexFiltering looks for a peak one isotope
    below) of isotope k of another channel, 0 <= k <= maxIsotope_(): a peak there belongs to the peptide, not to another envelope.
    Like MultiplexFiltering, which exempts only the isotopes of the pattern it detects, the higher isotopes of another channel do not
    count: e.g. the M-1 of a Lys8 channel lies 0.0063 Th (2+) from the light channel's M+7, which cannot explain a strong peak there.

    The decision is by theoretical position, while the peak read at the M-1 (lowerPeak_()) is the nearest actual centroid within the
    same window. Real isotope centroids are not spaced by the 13C mass difference (mainly their 15N and 34S variants pull M+3 to about
    +3.006-3.008 Da instead of +3.010), so another channel's isotope centroid can be read at the M-1 while its theoretical position
    lies just outside the window: in [][Lys4,Arg6][Lys8,Arg10], the medium M+3 at the heavy M-1 for peptides of about 1.7-2.1 kDa, where
    a heavy channel weaker than about the medium is then not quantifiable. Conversely, inside the window any peak at the M-1 is exempt,
    also an unrelated one. The lower isotopologue that a deuterium label (Lys4, Leu3, Dimethyl4/6/8, ICPL4/10) carries itself is not
    modelled either, as in MultiplexFiltering: the 2H3 form of a Lys4 (2H4) channel lies 1.0063/z Th below its M0, inside the window.
  */
  bool lowerIsChannel_(double lower_mz, int charge, const vector<double>& mzs, Size self) const
  {
    if (charge == 0) { return false; }
    const double spacing = Constants::C13C12_MASSDIFF_U / std::abs(charge);
    const double tolerance = 2.0 * detectionMZTolerance_(lower_mz);
    const double max_isotope = static_cast<double>(maxIsotope_());
    for (Size o = 0; o < mzs.size(); ++o)
    {
      if (o == self) { continue; }
      const double k = std::round((lower_mz - mzs[o]) / spacing);
      if (k >= 0.0 && k <= max_isotope && std::abs(lower_mz - mzs[o] - k * spacing) <= tolerance) { return true; }
    }
    return false;
  }

  /// The intensity of the nearest peak of @p spectrum to @p lower_mz (one isotope below a channel) within twice the detection's m/z tolerance, 0 if none
  double lowerPeak_(const MSSpectrum& spectrum, double lower_mz) const
  {
    const Int index = spectrum.findNearest(lower_mz, 2.0 * detectionMZTolerance_(lower_mz));
    return index >= 0 ? static_cast<double>(spectrum[index].getIntensity()) : 0.0;
  }

  /// The result of pairedIntensity_
  struct PairedIntensity_
  {
    double intensity = 0.0; ///< the channel's intensity, 0 if it shows no signal
    Size paired = 0;        ///< spectra in which at least one peak of the channel was paired with the reference
    Size pattern = 0;       ///< of them, the spectra that show the channel's isotope pattern
  };

  /**
    @brief The intensity of a channel from its peaks paired with those of a reference channel.

    @p channel and @p reference hold the intensities of the isotopic peaks of the two channels per spectrum ([isotope][spectrum],
    0 = no peak), @p lower the peak one isotope below the channel per spectrum (0 = none, or an isotope of another channel of the
    peptide: lowerIsChannel_()); spectra @p first to @p last are used. The intensity is @p reference_intensity scaled by the ratio
    of the summed channel and reference peaks of the same isotope and spectrum, i.e. the ratio is taken where both channels were
    measured, the way the detection corrects the channels of a multiplet. It is 0 (no signal) unless the first two isotopes of the
    channel (one if there is only one) occur together in at least two spectra.

    Also counted are the spectra in which at least one peak of the channel was paired with the reference, and of them the
    <em>pattern spectra</em>, which show the channel's isotope pattern: its first two isotopes (one if there is only one), and no
    peak one isotope lower that exceeds lower_isotope_fraction_ of its monoisotopic peak (else its peaks are the isotopes of
    another envelope, as in the detection). The callers accept the intensity only with at least min_pattern_spectra_ pattern
    spectra; with fewer, a channel with an intensity is not quantifiable, not absent. The intensity itself stays the one over all
    paired peaks: taken over the pattern spectra only, it was less accurate on real data (1:1 mixtures).
  */
  static PairedIntensity_ pairedIntensity_(const vector<vector<double>>& channel, const vector<vector<double>>& reference, const vector<double>& lower,
                                           Size first, Size last, double reference_intensity)
  {
    PairedIntensity_ result;
    const Size traces = std::min(channel.size(), reference.size());
    if (traces == 0 || channel[0].empty() || first > last || !(reference_intensity > 0.0)) { return result; }
    const Size required = std::min<Size>(traces, 2);
    double channel_sum = 0.0, reference_sum = 0.0;
    Size together = 0;
    for (Size i = first; i <= last; ++i)
    {
      Size found = 0;
      bool pair_seen = false;
      for (Size k = 0; k < traces; ++k)
      {
        if (!(channel[k][i] > 0.0)) { continue; }
        if (k < required) { ++found; }
        if (!(reference[k][i] > 0.0)) { continue; }
        channel_sum += channel[k][i];
        reference_sum += reference[k][i];
        pair_seen = true;
      }
      if (found == required) { ++together; }
      result.paired += pair_seen;
      if (found == required && pair_seen && !(lower[i] > lower_isotope_fraction_ * channel[0][i])) { ++result.pattern; }
    }
    if (together >= 2 && reference_sum > 0.0) { result.intensity = reference_intensity * channel_sum / reference_sum; }
    return result;
  }

  /**
    @brief Pearson correlation of the elution profiles of a channel and its reference

    @p channel and @p reference hold the intensities of their isotopic peaks per spectrum ([isotope][spectrum], 0 = no peak); a
    profile is the sum over the isotopes per spectrum. NaN with fewer than three spectra or a constant profile.
  */
  static double elutionCorrelation_(const vector<vector<double>>& channel, const vector<vector<double>>& reference)
  {
    const Size traces = std::min(channel.size(), reference.size());
    const Size n = traces == 0 ? 0 : std::min(channel[0].size(), reference[0].size());
    if (n < 3) { return std::numeric_limits<double>::quiet_NaN(); }
    vector<double> channel_profile(n, 0.0), reference_profile(n, 0.0);
    for (Size k = 0; k < traces; ++k)
    {
      for (Size i = 0; i < n; ++i)
      {
        channel_profile[i] += channel[k][i];
        reference_profile[i] += reference[k][i];
      }
    }
    return Math::pearsonCorrelationCoefficient(channel_profile.begin(), channel_profile.end(), reference_profile.begin(), reference_profile.end());
  }

  /// The result of requantifyChannel_
  struct Requantified_
  {
    double intensity = 0.0;  ///< > 0 measured, 0 absent, NaN not quantifiable
    Size paired = 0;         ///< spectra in which at least one peak of the channel was paired with the reference (whether accepted or not)
    bool no_pattern = false; ///< not quantifiable because its peaks do not show its isotope pattern in min_pattern_spectra_ spectra
    bool not_coeluting = false; ///< not quantifiable because its elution profile does not follow the reference's ('-requantify_min_correlation')
  };

  /**
    @brief The intensity of a channel that the multiplet detection did not find, measured at its expected position.

    Pairs the peaks of the first @p traces isotopic mass traces at @p mz (charge @p charge; one of the peptide's channels
    @p channel_mzs) with those of the same isotope of the reference channel at @p reference_mz, in the spectra of @p cv between
    @p rt_lo and @p rt_hi, within the detection's m/z tolerance (pairedIntensity_). (An integral of the weak channel on its own
    would miss the flanks and isotopes that fall below the intensity cutoff and be biased low.) A peak that a detected pattern
    reserved (@p blacklist) is not used for the channel if @p skip_reserved says so (it gets the blacklist spectrum and the index of
    the peak in it). The intensity is accepted if the channel's peaks show its isotope pattern paired with the reference in at
    least min_pattern_spectra_ spectra (the peak one isotope below it is read unless it is an isotope of another channel of
    @p channel_mzs, lowerIsChannel_()); otherwise the channel is not quantifiable (NaN, Requantified_::no_pattern), e.g. when
    its peaks are the isotopes of another envelope one isotope below it. Without an intensity, it is NaN (not quantifiable)
    instead of 0 if it would have one with the skipped peaks: they occupy its position, so its absence cannot be told.
  */
  Requantified_ requantifyChannel_(const MSExperiment& spectra, const MSExperiment& blacklist, const std::function<bool(const MSSpectrum&, Size)>& skip_reserved,
                                   const std::optional<double>& cv, double rt_lo, double rt_hi, double mz, int charge, const vector<double>& channel_mzs,
                                   Size traces, double reference_mz, double reference_intensity, bool check_coelution = false) const
  {
    Requantified_ result;
    if (charge == 0 || traces == 0 || !(reference_intensity > 0.0)) { return result; }
    const double spacing = Constants::C13C12_MASSDIFF_U / std::abs(charge);
    const double lower_mz = mz - spacing;
    const Size self = static_cast<Size>(std::find(channel_mzs.begin(), channel_mzs.end(), mz) - channel_mzs.begin());
    const bool lower_is_channel = lowerIsChannel_(lower_mz, charge, channel_mzs, self);
    vector<vector<double>> channel(traces), occupied(traces), reference(traces); // occupied: with the skipped peaks
    vector<double> lower; // per spectrum: the peak one isotope below the channel where its monoisotopic peak is found
    bool skipped = false;
    for (auto it = spectra.RTBegin(rt_lo); it != spectra.RTEnd(rt_hi); ++it)
    {
      if (compensationVoltage_(*it) != cv) { continue; }
      const MSSpectrum* black = blacklistSpectrum_(blacklist, *it, cv);
      double lower_intensity = 0.0;
      for (Size k = 0; k < traces; ++k)
      {
        const double shift = static_cast<double>(k) * spacing;
        double channel_intensity = 0.0, occupied_intensity = 0.0, reference_intensity_k = 0.0;
        const Int index = it->findNearest(mz + shift, detectionMZTolerance_(mz + shift));
        if (index >= 0)
        {
          occupied_intensity = (*it)[index].getIntensity();
          const Int reserved = black != nullptr ? black->findNearest((*it)[index].getMZ(), 1e-9) : -1;
          if (reserved >= 0 && skip_reserved(*black, static_cast<Size>(reserved))) { skipped = true; }
          else { channel_intensity = occupied_intensity; }
          const Int ref = it->findNearest(reference_mz + shift, detectionMZTolerance_(reference_mz + shift));
          if (ref >= 0) { reference_intensity_k = (*it)[ref].getIntensity(); }
          if (k == 0 && !lower_is_channel) { lower_intensity = lowerPeak_(*it, lower_mz); }
        }
        channel[k].push_back(channel_intensity);
        occupied[k].push_back(occupied_intensity);
        reference[k].push_back(reference_intensity_k);
      }
      lower.push_back(lower_intensity);
    }
    if (channel[0].empty()) { return result; }
    const Size last = channel[0].size() - 1;
    const PairedIntensity_ paired = pairedIntensity_(channel, reference, lower, 0, last, reference_intensity);
    result.paired = paired.paired;
    if (paired.intensity > 0.0 && paired.pattern >= min_pattern_spectra_)
    {
      // the channels of one peptide co-elute: a weak signal whose elution profile does not follow the reference's is another species
      if (check_coelution && requantify_min_correlation_ > -1.0 && paired.intensity < requantify_correlation_ratio_ * reference_intensity
          && !(elutionCorrelation_(channel, reference) >= requantify_min_correlation_))
      {
        result.intensity = std::numeric_limits<double>::quiet_NaN();
        result.not_coeluting = true;
      }
      else { result.intensity = paired.intensity; }
    }
    else if (paired.intensity > 0.0)
    {
      result.intensity = std::numeric_limits<double>::quiet_NaN();
      result.no_pattern = true;
    }
    else if (skipped && pairedIntensity_(occupied, reference, lower, 0, last, reference_intensity).intensity > 0.0)
    {
      result.intensity = std::numeric_limits<double>::quiet_NaN();
    }
    return result;
  }

  /**
    @brief Re-quantify the channels that completed multiplets lack (dummy features with intensity 0).

    A partial multiplet is completed by the resolver with a dummy feature per missing channel: intensity NaN if another
    detected pattern overlaps with its position (not quantifiable), otherwise 0. Such a 0 only says that the detection did
    not find the channel, which can be too weak for its isotope pattern, correlation or intensity filters (e.g. a 1:10
    mixture). It is replaced by the signal at the channel's expected position over the elution region of the multiplet,
    relative to the most intense detected channel and with as many mass traces as the detected channels have
    (requantifyChannel_), so that a weak channel still yields a ratio and only a really absent one is reported as 0 (like
    MaxQuant's 're-quantify'). Of the reserved peaks, only those of no reported feature are used if the blacklist records
    @p owners (mayUseReserved_(): the missing channel was not detected, so the multiplet's own reserved peaks at its position
    are isotopes of its other channels), none without. A channel that has no value only because such peaks occupy its
    position is not quantifiable (NaN), not absent; so is a channel whose peaks do not show its isotope pattern paired with the
    reference in at least min_pattern_spectra_ spectra (e.g. the isotopes of another envelope one isotope below it). The peak one
    isotope below a channel is not taken for another envelope if it is an isotope of another channel of the multiplet
    (lowerIsChannel_(), with the m/z of all its channels).
  */
  void requantifyMissingChannels_(ConsensusMap& resolved, const FeatureMap& channels, const MSExperiment& centroided,
                                  const MSExperiment& blacklist, bool owners) const
  {
    std::unordered_map<UInt64, const Feature*> channel_of;
    for (const Feature& channel : channels) { channel_of.emplace(channel.getUniqueId(), &channel); }

    Size requantified = 0, absent = 0, occupied = 0, no_pattern = 0, not_coeluting = 0, reserved_used = 0;
    for (auto& cf : resolved)
    {
      // the elution region and the number of mass traces of the detected channels, and the most intense one
      double rt_lo = std::numeric_limits<double>::max(), rt_hi = std::numeric_limits<double>::lowest();
      Size traces = 0;
      bool has_dummy = false;
      const FeatureHandle* reference = nullptr;
      for (const FeatureHandle& handle : cf.getFeatures())
      {
        const auto it = channel_of.find(handle.getUniqueId());
        if (it == channel_of.end())
        {
          has_dummy = has_dummy || handle.getIntensity() == 0.0;
          continue;
        }
        const DBoundingBox<2> box = it->second->getConvexHull().getBoundingBox();
        rt_lo = std::min(rt_lo, box.minX());
        rt_hi = std::max(rt_hi, box.maxX());
        traces = std::max(traces, it->second->getConvexHulls().size());
        if (reference == nullptr || handle.getIntensity() > reference->getIntensity()) { reference = &handle; }
      }
      if (!has_dummy || traces == 0 || reference == nullptr) { continue; }
      const double reference_mz = reference->getMZ(), reference_intensity = reference->getIntensity();

      const std::optional<double> cv = compensationVoltage_(cf);
      auto skip_reserved = [&](const MSSpectrum& black, Size index)
      {
        if (!owners || !mayUseReserved_(MultiplexResolverAlgorithm::blacklistOwner(black, index), UniqueIdInterface::INVALID)) { return true; }
        ++reserved_used;
        return false;
      };
      vector<double> channel_mzs; // all channels of the multiplet, the missing ones included
      for (const FeatureHandle& handle : cf.getFeatures()) { channel_mzs.push_back(handle.getMZ()); }
      ConsensusFeature::HandleSetType handles;
      IntList requantified_channels; // channel indices (the columns of this run's map) that got a value
      for (FeatureHandle handle : cf.getFeatures())
      {
        if (handle.getIntensity() == 0.0 && !channel_of.contains(handle.getUniqueId()))
        {
          const Requantified_ result = requantifyChannel_(centroided, blacklist, skip_reserved, cv, rt_lo, rt_hi, handle.getMZ(), handle.getCharge(),
                                                          channel_mzs, traces, reference_mz, reference_intensity, true);
          const double intensity = result.intensity;
          if (intensity > 0.0 || std::isnan(intensity)) { handle.setIntensity(intensity); }
          requantified += intensity > 0.0;
          if (intensity > 0.0) { requantified_channels.push_back(static_cast<Int>(handle.getMapIndex())); }
          occupied += std::isnan(intensity) && !result.no_pattern && !result.not_coeluting;
          no_pattern += result.no_pattern;
          not_coeluting += result.not_coeluting;
          absent += intensity == 0.0;
        }
        handles.insert(handle);
      }
      cf.setFeatures(std::move(handles));
      // the handles lose their run in linking, the identifications keep it (map_index): see markRequantified_()
      if (!requantified_channels.empty())
      {
        for (PeptideIdentification& id : cf.getPeptideIdentifications()) { id.setMetaValue(requantified_channels_meta_, requantified_channels); }
      }
    }
    if (requantified + absent + occupied + no_pattern + not_coeluting > 0)
    {
      OPENMS_LOG_INFO << "Re-quantified " << requantified << " channel(s) the detection did not find; " << absent
                      << " show no signal at their expected position and are reported as absent, " << occupied
                      << " lie on reserved peaks they may not use (e.g. isotopes of the multiplet's other channels) and are not quantifiable." << endl;
    }
    if (no_pattern > 0)
    {
      OPENMS_LOG_INFO << no_pattern << " re-quantified channel(s) are not quantifiable: their peaks do not show the channel's isotope pattern (both "
                      << "of its first two isotopes, and no peak one isotope lower that exceeds 30 % of its monoisotopic peak) in at least "
                      << min_pattern_spectra_ << " spectra paired with the reference." << endl;
    }
    if (not_coeluting > 0)
    {
      OPENMS_LOG_INFO << not_coeluting << " re-quantified channel(s) weaker than " << requantify_correlation_ratio_ << " times the reference are not "
                      << "quantifiable: their elution profile correlates with the reference's by less than " << requantify_min_correlation_
                      << " ('-requantify_min_correlation')." << endl;
    }
    if (reserved_used > 0)
    {
      OPENMS_LOG_INFO << "The re-quantification used " << reserved_used << " reserved peak(s) of no reported feature (patterns the detection dropped)." << endl;
    }
  }

  /// Isotopic peaks per channel that the targeted quantification uses (the minimum of 'algorithm:isotopes_per_peptide' by default)
  static constexpr Size targeted_isotopes_ = 3;
  /// Spectra in which the identified channel needs its first two isotopes for a targeted feature
  static constexpr Size min_targeted_reference_spectra_ = 3;
  /// A channel's peak is left out (the identified channel's: its intensity is not reported) if isotopes of other channels are expected to make up at least this fraction of it
  static constexpr double targeted_other_isotope_fraction_ = 0.2;

  /// An unassigned identification that the targeted quantification may measure
  struct TargetedCandidate_
  {
    Size index = 0;             ///< position among the unassigned identifications
    std::string peptide;        ///< MultiplexResolverAlgorithm::peptideKey() of its best hit
    vector<double> mzs;         ///< expected m/z of the peptide per channel
    Size own = 0;               ///< channel of the identification
    int charge = 0;
    double rt = 0.0;
    std::optional<double> cv;
    vector<char> on_other_multiplet; ///< per channel: lies on a channel of another peptide's multiplet (not quantifiable)
  };

  /// One elution peak measured by the targeted quantification
  struct TargetedPeak_
  {
    double rt_lo = 0.0, rt_hi = 0.0, rt = 0.0; ///< boundaries, and the RT of the identified channel
    vector<double> intensities;                 ///< per channel: > 0 measured, 0 absent, NaN not quantifiable
    vector<UInt64> owners;                      ///< features whose reserved peaks entered its reported intensities
  };

  enum class TargetedOutcome_ { MEASURED, NO_SIGNAL, ISOTOPE_PATTERN };

  /**
    @brief Measure the peptide of @p candidate (identified as @p hit) in the MS1 peaks around its identification.

    See quantifyTargeted_(). @p occupied tells whether a reserved peak (blacklist spectrum, index) belongs to another peptide.
  */
  TargetedOutcome_ measureTargeted_(const TargetedCandidate_& candidate, const PeptideHit& hit, const MSExperiment& centroided,
                                    const MSExperiment& blacklist, const std::function<bool(const MSSpectrum&, Size)>& occupied,
                                    ElutionPeakDetection& elution, const PeakIntegrator& integrator, TargetedPeak_& peak) const
  {
    const double rt_typical = getParam_().getValue("algorithm:rt_typical");
    const double rt_tolerance = getParam_().getValue("id_mapping:rt_tolerance");
    const Size n_channels = candidate.mzs.size(), own = candidate.own;
    const double spacing = Constants::C13C12_MASSDIFF_U / std::abs(candidate.charge);

    // A strong peak one isotope below a channel is the envelope of another species (unless it is an isotope of another channel,
    // lowerIsChannel_()), whose isotopes the channel's peaks would be.
    vector<char> lower_is_channel(n_channels);
    for (Size c = 0; c < n_channels; ++c) { lower_is_channel[c] = lowerIsChannel_(candidate.mzs[c] - spacing, candidate.charge, candidate.mzs, c); }

    // A peak of channel c can hold an isotope (k >= 1, up to the isotopes the detection assigns to a channel) of a lower channel o
    // (see quantifyTargeted_()): it does if the peak is also the nearest one to the isotope's m/z within twice the detection's m/z
    // tolerance (the window in which MultiplexFiltering looks for peaks of another envelope), i.e. the two are one centroid; a
    // resolved isotope has a centroid of its own, nearer to its m/z. Its part is o's monoisotopic peak, without the isotopes of the
    // channels below o (hence the correction in ascending m/z), times I_k/I_0 of the peptide's isotope distribution.
    const Size max_isotope = maxIsotope_();
    const IsotopeDistribution theoretical = hit.getSequence().getFormula().getIsotopeDistribution(
      CoarseIsotopePatternGenerator(std::max(targeted_isotopes_, max_isotope + 1)));
    vector<Size> ascending_mz(n_channels);
    std::iota(ascending_mz.begin(), ascending_mz.end(), Size(0));
    std::sort(ascending_mz.begin(), ascending_mz.end(), [&](Size a, Size b) { return candidate.mzs[a] < candidate.mzs[b]; });
    // the expected part of the lower channels' isotopes in peak @p index of @p spectrum, found for channel @p c
    auto isotopes_of_others = [&](const MSSpectrum& spectrum, Size c, Int index, const vector<double>& monoisotopic)
    {
      double expected = 0.0;
      for (Size o = 0; o < n_channels; ++o)
      {
        const double k = std::round((spectrum[index].getMZ() - candidate.mzs[o]) / spacing);
        if (o == c || !(monoisotopic[o] > 0.0) || k < 1.0 || k > static_cast<double>(max_isotope) || static_cast<Size>(k) >= theoretical.size())
        {
          continue;
        }
        const double isotope_mz = candidate.mzs[o] + k * spacing;
        if (spectrum.findNearest(isotope_mz, 2.0 * detectionMZTolerance_(isotope_mz)) != index) { continue; }
        expected += monoisotopic[o] * theoretical[static_cast<Size>(k)].getIntensity() / theoretical[0].getIntensity();
      }
      return expected;
    };

    // per channel, isotope and spectrum of the window: the nearest peak (0 if none, or if it belongs to another peptide or channel),
    // without the expected part of the other channels' isotopes (at least 0); a peak of another channel than the identified one is
    // left out if that part is at least targeted_other_isotope_fraction_ of it, the identified channel's peaks (the reference) are
    // kept, and also kept as measured (own_measured)
    vector<double> rts;
    vector<vector<vector<double>>> intensity(n_channels, vector<vector<double>>(targeted_isotopes_));
    vector<vector<double>> own_measured(targeted_isotopes_); // the identified channel's peaks without the subtraction
    vector<vector<char>> reserved(n_channels); // a peak of an isotope read (M+0 to M+2) belongs to another peptide or channel (left out)
    struct ReadOwner { Size spectrum, channel, isotope; UInt64 owner; };
    vector<ReadOwner> read_owners; // the reserved peaks read
    vector<vector<double>> lower_peaks(n_channels); // per channel: the peak one isotope below it, unless that is an isotope of another channel
    vector<char> own_mostly_other; // a peak of the identified channel's first two isotopes is mostly another channel's isotope
    vector<double> monoisotopic(n_channels); // per channel: its peak at its m/z in the spectrum without the lower channels' isotopes
    for (auto it = centroided.RTBegin(candidate.rt - rt_typical); it != centroided.RTEnd(candidate.rt + rt_typical); ++it)
    {
      if (compensationVoltage_(*it) != candidate.cv) { continue; }
      const MSSpectrum* black = blacklistSpectrum_(blacklist, *it, candidate.cv);
      rts.push_back(it->getRT());
      std::fill(monoisotopic.begin(), monoisotopic.end(), 0.0);
      for (Size c : ascending_mz)
      {
        const Int index = it->findNearest(candidate.mzs[c], detectionMZTolerance_(candidate.mzs[c]));
        if (index >= 0) { monoisotopic[c] = std::max(0.0, (*it)[index].getIntensity() - isotopes_of_others(*it, c, index, monoisotopic)); }
      }
      own_mostly_other.push_back(false);
      for (Size c = 0; c < n_channels; ++c)
      {
        reserved[c].push_back(false);
        for (Size k = 0; k < targeted_isotopes_; ++k)
        {
          const double mz = candidate.mzs[c] + static_cast<double>(k) * spacing;
          const Int index = it->findNearest(mz, detectionMZTolerance_(mz));
          const Int reserved_index = (index >= 0 && black != nullptr) ? black->findNearest((*it)[index].getMZ(), 1e-9) : -1;
          double value = index >= 0 ? (*it)[index].getIntensity() : 0.0;
          bool other = index >= 0
                       && ((c != own && candidate.on_other_multiplet[c])
                           || (reserved_index >= 0 && occupied(*black, static_cast<Size>(reserved_index))));
          if (index >= 0 && !other)
          {
            const double expected = isotopes_of_others(*it, c, index, monoisotopic);
            const bool mostly_other = expected > 0.0 && expected >= targeted_other_isotope_fraction_ * value;
            value = std::max(0.0, value - expected);
            if (c != own) { other = mostly_other; }
            else if (mostly_other && k < 2) { own_mostly_other.back() = true; }
          }
          intensity[c][k].push_back(index >= 0 && !other ? value : 0.0);
          if (c == own) { own_measured[k].push_back(index >= 0 && !other ? (*it)[index].getIntensity() : 0.0); }
          if (other) { reserved[c].back() = true; }
          if (reserved_index >= 0 && !other)
          {
            read_owners.push_back({rts.size() - 1, c, k, MultiplexResolverAlgorithm::blacklistOwner(*black, static_cast<Size>(reserved_index))});
          }
        }
      }
      for (Size c = 0; c < n_channels; ++c) { lower_peaks[c].push_back(lower_is_channel[c] ? 0.0 : lowerPeak_(*it, candidate.mzs[c] - spacing)); }
    }
    if (rts.size() < min_targeted_reference_spectra_) { return TargetedOutcome_::NO_SIGNAL; }

    // The elution peaks of the sum over the channels' first two isotopes, smoothed and split at minima as ElutionPeakDetection does
    // (its detectPeaks() would also build sub-traces, which fails for one of zeros only), each trimmed to where the smoothed sum is
    // positive around its apex. The peptide's is the one at the identification.
    vector<Peak2D> points(rts.size());
    for (Size i = 0; i < rts.size(); ++i)
    {
      double sum = 0.0;
      for (Size c = 0; c < n_channels; ++c) { sum += intensity[c][0][i] + intensity[c][1][i]; }
      points[i].setRT(rts[i]);
      points[i].setMZ(candidate.mzs[own]);
      points[i].setIntensity(sum);
    }
    MassTrace trace(points);
    const double chrom_fwhm = elution.getParameters().getValue("chrom_fwhm");
    const Size window = static_cast<Size>(std::ceil(chrom_fwhm / trace.getAverageMS1CycleTime()));
    elution.smoothData(trace, static_cast<int>(window));
    vector<Size> maxima, minima;
    elution.findLocalExtrema(trace, window / 2, maxima, minima);
    const vector<double>& smoothed = trace.getSmoothedIntensities();
    minima.push_back(rts.size() - 1); // the end of the last part
    Size first = 0, last = 0;
    double best_distance = std::numeric_limits<double>::infinity();
    for (Size part = 0, begin = 0; part < minima.size(); begin = minima[part] + 1, ++part)
    {
      if (begin > minima[part]) { continue; }
      const Size top = static_cast<Size>(std::max_element(smoothed.begin() + begin, smoothed.begin() + minima[part] + 1) - smoothed.begin());
      if (!(smoothed[top] > 0.0)) { continue; }
      Size lo = top, hi = top;
      while (lo > begin && smoothed[lo - 1] > 0.0) { --lo; }
      while (hi < minima[part] && smoothed[hi + 1] > 0.0) { ++hi; }
      const double distance = std::max({0.0, rts[lo] - candidate.rt, candidate.rt - rts[hi]});
      if (distance <= rt_tolerance && distance < best_distance)
      {
        best_distance = distance;
        first = lo;
        last = hi;
      }
    }
    if (std::isinf(best_distance)) { return TargetedOutcome_::NO_SIGNAL; }

    // the identified channel is the reference, measured where its first two isotopes are seen without another envelope below (a
    // strong peak one isotope below). If the subtraction of the other channels' isotopes leaves it fewer than
    // min_targeted_reference_spectra_ such spectra although its peaks as measured have them, the subtracted part cannot be right,
    // since the channel is identified there (e.g. an unrelated peak at a lower channel's m/z was taken for that channel): its peaks
    // are then used as measured, and its own intensity is not quantifiable.
    vector<vector<double>> reference;
    Size clean = 0, mostly_own = 0; // mostly_own: of them, spectra in which its first two isotopes are not mostly another channel's isotope
    auto take_reference = [&](const vector<vector<double>>& peaks)
    {
      reference.assign(targeted_isotopes_, vector<double>(rts.size(), 0.0));
      clean = 0;
      mostly_own = 0;
      for (Size i = first; i <= last; ++i)
      {
        if (!(peaks[0][i] > 0.0) || !(peaks[1][i] > 0.0) || lower_peaks[own][i] > lower_isotope_fraction_ * peaks[0][i]) { continue; }
        ++clean;
        mostly_own += !own_mostly_other[i];
        for (Size k = 0; k < targeted_isotopes_; ++k) { reference[k][i] = peaks[k][i]; }
      }
    };
    take_reference(intensity[own]);
    const bool own_as_measured = clean < min_targeted_reference_spectra_;
    if (own_as_measured) { take_reference(own_measured); }
    if (clean < min_targeted_reference_spectra_) { return TargetedOutcome_::NO_SIGNAL; }
    const vector<vector<double>>& own_peaks = own_as_measured ? own_measured : intensity[own];

    // its isotope pattern has to fit the peptide (the detection's averagine filter, on the identified sequence)
    vector<double> model(targeted_isotopes_, 0.0), measured(targeted_isotopes_, 0.0);
    for (Size k = 0; k < targeted_isotopes_; ++k)
    {
      if (k < theoretical.size()) { model[k] = theoretical[k].getIntensity(); }
      measured[k] = std::accumulate(reference[k].begin(), reference[k].end(), 0.0);
    }
    const double averagine_similarity = getParam_().getValue("algorithm:averagine_similarity");
    if (!(Math::pearsonCorrelationCoefficient(model.begin(), model.end(), measured.begin(), measured.end()) >= averagine_similarity))
    {
      return TargetedOutcome_::ISOTOPE_PATTERN;
    }

    // its intensity and RT as the detection determines them: the trapezoid areas of its isotopic traces, and the intensity-weighted
    // mean RT of its monoisotopic peaks
    double reference_intensity = 0.0, weighted_rt = 0.0;
    for (Size i = first; i <= last; ++i) { weighted_rt += reference[0][i] * rts[i]; }
    peak.rt = weighted_rt / std::accumulate(reference[0].begin() + first, reference[0].begin() + last + 1, 0.0);
    for (Size k = 0; k < targeted_isotopes_; ++k)
    {
      MSChromatogram chromatogram;
      for (Size i = first; i <= last; ++i)
      {
        if (reference[k][i] > 0.0) { chromatogram.push_back(ChromatogramPeak(rts[i], reference[k][i])); }
      }
      if (chromatogram.size() > 1) { reference_intensity += integrator.integratePeak(chromatogram, chromatogram.front().getRT(), chromatogram.back().getRT()).area; }
      else if (chromatogram.size() == 1) { reference_intensity += chromatogram[0].getIntensity(); }
    }

    // every other channel from its peaks paired with the reference, if they show its isotope pattern and it co-elutes with the reference
    // (compared with all peaks of the identified channel: its peaks in the spectra left out above still show where it elutes)
    const double peptide_similarity = getParam_().getValue("algorithm:peptide_similarity");
    auto profile = [&](const vector<vector<double>>& peaks)
    {
      vector<double> summed(last - first + 1, 0.0);
      for (Size k = 0; k < targeted_isotopes_; ++k)
      {
        for (Size i = first; i <= last; ++i) { summed[i - first] += peaks[k][i]; }
      }
      return summed;
    };
    const vector<double> reference_profile = profile(own_peaks);
    peak.intensities.assign(n_channels, 0.0);
    // not quantifiable itself if its peaks are mostly another channel's isotope (as another channel's would be left out)
    peak.intensities[own] = !own_as_measured && mostly_own >= min_targeted_reference_spectra_ ? reference_intensity : std::numeric_limits<double>::quiet_NaN();
    for (Size c = 0; c < n_channels; ++c)
    {
      if (c == own) { continue; }
      if (candidate.on_other_multiplet[c])
      {
        peak.intensities[c] = std::numeric_limits<double>::quiet_NaN();
        continue;
      }
      const PairedIntensity_ paired = pairedIntensity_(intensity[c], reference, lower_peaks[c], first, last, reference_intensity);
      const vector<double> channel_profile = profile(intensity[c]);
      const double correlation = Math::pearsonCorrelationCoefficient(channel_profile.begin(), channel_profile.end(), reference_profile.begin(), reference_profile.end());
      if (paired.intensity > 0.0 && paired.pattern >= min_pattern_spectra_ && correlation >= peptide_similarity)
      {
        peak.intensities[c] = paired.intensity;
        continue;
      }
      // not quantifiable if something is there: a signal that does not co-elute or does not show the channel's isotope pattern, a
      // peak left out above (of another peptide or channel) at any isotope read (M+0 to M+2), or any other peak of the channel (its
      // absence cannot be told); absent (0) only if none of these positions holds a peak, its own or another's, within the elution peak
      const bool interfered = paired.intensity > 0.0 || std::any_of(reserved[c].begin() + first, reserved[c].begin() + last + 1, [](char r) { return r != 0; });
      const bool seen = std::any_of(intensity[c].begin(), intensity[c].end(), [&](const vector<double>& isotope)
                                    { return std::any_of(isotope.begin() + first, isotope.begin() + last + 1, [](double v) { return v > 0.0; }); });
      peak.intensities[c] = interfered || seen ? std::numeric_limits<double>::quiet_NaN() : 0.0;
    }
    peak.rt_lo = rts[first];
    peak.rt_hi = rts[last];
    // the owners of the reserved peaks that entered a reported intensity (the reference, or a channel measured from paired peaks)
    for (const ReadOwner& read : read_owners)
    {
      const Size i = read.spectrum;
      if (i >= first && i <= last && read.owner != UniqueIdInterface::INVALID && reference[read.isotope][i] > 0.0
          && intensity[read.channel][read.isotope][i] > 0.0 && (read.channel == own || peak.intensities[read.channel] > 0.0))
      {
        peak.owners.push_back(read.owner);
      }
    }
    return TargetedOutcome_::MEASURED;
  }

  /**
    @brief Quantify identified peptides that have no multiplet in the run from the MS1 peaks where they are expected ('-targeted_quantification').

    The candidates are the identifications that are not attached to a resolved multiplet: target hits whose labels fit a pattern
    (MultiplexResolverAlgorithm::expectedChannelMZs()), of a peptide (sequence without labels, charge, CV) that no resolved
    multiplet quantifies within 'algorithm:rt_typical' of its elution region - the abundances of a peptide add up over its
    multiplets, so it must not be measured twice. Nor is an identification whose own channel lies at a quantified channel of a
    resolved multiplet of another peptide (same charge and CV), within 'algorithm:rt_typical' of its elution region: its signal is
    that multiplet's (e.g. its tail beyond the reserved peaks), whose identification the resolver preferred. Another channel there
    is not quantifiable.

    Around each identification (+/- 'algorithm:rt_typical'), the first isotopic peaks of every channel are taken from the spectra
    the detection worked on (nearest peak within 'algorithm:mz_tolerance', same CV). Left out are the peaks the detection reserved
    for a multiplet of another identified peptide (it is quantified there) and the peaks of an unidentified multiplet (of
    @p multiplets, the detected ones) that cannot be the peptide (other charge or CV, or a channel not at one of the peptide's).
    Unlike in the re-quantification, the peaks of other unidentified multiplets are used, since they can be the peptide itself,
    detected without its identification mapping onto it.

    A channel's peak can also hold an isotope (M+1 up to M+n, n + 1 the maximum of 'algorithm:isotopes_per_peptide') of a lower
    channel, e.g. a Dimethyl4 or Lys4 channel the light channel's M+4 (0.0117 Da away), a Dimethyl8 channel the Dimethyl4 channel's
    M+4 (0.0059 Da). The isotope counts as part of the peak if the peak is also the nearest one to the isotope's m/z within twice
    'algorithm:mz_tolerance' (the window in which the detection looks for peaks of another envelope): the two are one centroid,
    unresolved, while a resolved isotope has a centroid of its own. Its part is the lower channel's monoisotopic peak in the
    spectrum times I_k/I_0 of the identified sequence's isotope distribution; that monoisotopic peak is itself taken without the
    isotopes of the channels below it (the channels are corrected in ascending m/z, so that e.g. the light M+4 in the Dimethyl4
    peak is not taken for Dimethyl4 when its M+4 is taken off the Dimethyl8 peak). The part is subtracted from the peak; a peak of
    another channel than the identified one is left out instead if the part is at least targeted_other_isotope_fraction_ of it
    (the estimate's error would be a substantial part of what is left, e.g. a Lys4 channel on the M+4 of a much stronger light
    channel). The elution peak is taken from the sum of the channels' first two isotopes (ElutionPeakDetection: smoothed and split
    at minima): the one containing the identification, else the nearest within 'id_mapping:rt_tolerance'.

    Within its boundaries, the identified channel is the reference. It needs its first two isotopes in at least
    min_targeted_reference_spectra_ spectra (not counting those with a strong peak one isotope below, i.e. another envelope, as
    in the detection, unless that peak is an isotope of another channel: lowerIsChannel_()) and an isotope pattern that fits the
    identified sequence (Pearson >= 'algorithm:averagine_similarity');
    its intensity is the trapezoid area of its isotopic traces, as the detection integrates them. Its peaks keep their place as
    the reference also where another channel's isotope makes up at least targeted_other_isotope_fraction_ of them (only the
    subtracted part's error is then not small; the other channels' paired intensities scale with the reference and do not depend
    on a constant error of it), but its own intensity is then reported not quantifiable, unless min_targeted_reference_spectra_ of
    its spectra are free of that, so that the same data give the same values whichever channel was identified. If the subtraction
    leaves it fewer than min_targeted_reference_spectra_ spectra with its first two isotopes although its peaks as measured have
    them, the subtracted part cannot be right (the channel is identified there; e.g. an unrelated peak at a lower channel's m/z was
    taken for that channel): its peaks are then used as measured, and its own intensity is not quantifiable. Every other channel
    is measured against it from paired peaks (pairedIntensity_), and reported if its peaks show its isotope pattern paired with the
    reference in at least min_pattern_spectra_ spectra (first two isotopes, no strong peak one isotope below, as for the
    reference) and it co-elutes (Pearson of the isotope-summed profiles >= 'algorithm:peptide_similarity'). Otherwise it is
    absent (0) only if none of its isotopic peaks read (M+0 to M+2) is seen within the boundaries; it is not quantifiable (left
    out) if its peaks there do not form its isotope pattern, or if something else, or a peak left out above, is at any of them.

    Targeted peaks of different peptides that read the same peaks (the identified channel of one at a measured channel or at the
    identified channel of the other - the reference, also if its own intensity is not quantifiable -, same charge and CV,
    overlapping boundaries; e.g. isomers, or a wrong identification next to a right one) are one signal:
    only the one with the best identification is kept, as IDConflictResolverAlgorithm keeps one identification per multiplet.

    Each elution peak gives one consensus feature (RT: the intensity-weighted mean of the identified channel's monoisotopic peaks;
    m/z: the lightest channel's), onto which every identification of the peptide within the peak is moved, with the map index of
    its own channel. The feature and these identifications are marked with the meta value 'quantification_method' = 'targeted'
    (only the identifications' survives linking). An unidentified multiplet kept for match between runs whose reserved peaks
    entered a feature's intensities is that feature's signal and is removed, so that linking cannot quantify it a second time.
  */
  void quantifyTargeted_(ConsensusMap& resolved, const ConsensusMap& multiplets, const FeatureMap& channels, const MultiplexResolverAlgorithm& resolver,
                         const MSExperiment& centroided, const MSExperiment& blacklist, const MultiplexResolverAlgorithm::OwnerPeptides* peptides) const
  {
    const double rt_typical = getParam_().getValue("algorithm:rt_typical");
    const double rt_tolerance = getParam_().getValue("id_mapping:rt_tolerance");
    std::unordered_map<UInt64, const Feature*> channel_of;
    for (const Feature& channel : channels) { channel_of.emplace(channel.getUniqueId(), &channel); }

    // elution regions of the peptides (and CVs) that resolved multiplets quantify, and the channels they quantify there
    using Key = std::pair<std::string, std::optional<double>>;
    std::map<Key, vector<std::pair<double, double>>> quantified;
    struct QuantifiedChannel
    {
      double mz;
      Key key;
      int charge;
      double rt_lo, rt_hi;
    };
    vector<QuantifiedChannel> quantified_channels; // sorted by m/z
    for (const auto& cf : resolved)
    {
      if (cf.getPeptideIdentifications().empty() || cf.getPeptideIdentifications()[0].getHits().empty()) { continue; }
      double rt_lo = cf.getRT(), rt_hi = cf.getRT();
      for (const FeatureHandle& handle : cf.getFeatures())
      {
        if (const auto it = channel_of.find(handle.getUniqueId()); it != channel_of.end())
        {
          const DBoundingBox<2> box = it->second->getConvexHull().getBoundingBox();
          rt_lo = std::min(rt_lo, box.minX());
          rt_hi = std::max(rt_hi, box.maxX());
        }
      }
      Key key{resolver.peptideKey(cf.getPeptideIdentifications()[0].getHits()[0], cf.getCharge()), compensationVoltage_(cf)};
      for (const FeatureHandle& handle : cf.getFeatures())
      {
        if (handle.getIntensity() > 0.0) { quantified_channels.push_back({handle.getMZ(), key, cf.getCharge(), rt_lo, rt_hi}); }
      }
      quantified[std::move(key)].emplace_back(rt_lo, rt_hi);
    }
    std::sort(quantified_channels.begin(), quantified_channels.end(), [](const auto& a, const auto& b) { return a.mz < b.mz; });

    PeptideIdentificationList& unassigned = resolved.getUnassignedPeptideIdentifications();
    Size decoys = 0, no_pattern = 0, already = 0, other_multiplet = 0;
    std::map<Key, vector<TargetedCandidate_>> groups;
    for (Size i = 0; i < unassigned.size(); ++i)
    {
      const PeptideIdentification& id = unassigned[i];
      if (id.getHits().empty() || id.getHits()[0].getCharge() == 0) { ++no_pattern; continue; }
      const PeptideHit& hit = id.getHits()[0];
      if (hit.isDecoy()) { ++decoys; continue; }
      TargetedCandidate_ candidate;
      candidate.mzs = resolver.expectedChannelMZs(hit.getSequence(), hit.getCharge(), candidate.own);
      if (candidate.mzs.empty()) { ++no_pattern; continue; }
      candidate.index = i;
      candidate.charge = hit.getCharge();
      candidate.rt = id.getRT();
      candidate.cv = compensationVoltage_(id);
      candidate.peptide = resolver.peptideKey(hit, candidate.charge);
      Key key{candidate.peptide, candidate.cv};
      auto near = [&](double lo, double hi) { return candidate.rt >= lo - rt_typical && candidate.rt <= hi + rt_typical; };
      if (const auto it = quantified.find(key); it != quantified.end()
          && std::any_of(it->second.begin(), it->second.end(), [&](const auto& region) { return near(region.first, region.second); }))
      {
        ++already;
        continue;
      }
      // a channel on a channel of another peptide's multiplet, near its elution region, holds that multiplet's signal: the
      // identification is not measured if it is the identified channel, else that channel is not quantifiable
      auto on_other_multiplet = [&](double mz)
      {
        const double tolerance = detectionMZTolerance_(mz);
        for (auto it = std::lower_bound(quantified_channels.begin(), quantified_channels.end(), mz - tolerance,
                                        [](const QuantifiedChannel& q, double value) { return q.mz < value; });
             it != quantified_channels.end() && it->mz <= mz + tolerance; ++it)
        {
          if (it->charge == candidate.charge && it->key.second == candidate.cv && it->key.first != candidate.peptide && near(it->rt_lo, it->rt_hi)) { return true; }
        }
        return false;
      };
      if (on_other_multiplet(candidate.mzs[candidate.own]))
      {
        ++other_multiplet;
        continue;
      }
      candidate.on_other_multiplet.assign(candidate.mzs.size(), false);
      for (Size c = 0; c < candidate.mzs.size(); ++c) { candidate.on_other_multiplet[c] = c != candidate.own && on_other_multiplet(candidate.mzs[c]); }
      groups[std::move(key)].push_back(std::move(candidate));
    }

    ElutionPeakDetection elution;
    {
      Param elution_param = elution.getDefaults();
      elution_param.setValue("width_filtering", "off");
      elution.setParameters(elution_param);
    }
    PeakIntegrator integrator;
    {
      Param integrator_param = integrator.getDefaults();
      integrator_param.setValue("integration_type", "trapezoid");
      integrator.setParameters(integrator_param);
    }

    // the elution peaks of every peptide, with their identifications (positions in the group)
    struct Measured
    {
      const vector<TargetedCandidate_>* group;
      TargetedPeak_ peak;
      vector<Size> members;
    };
    vector<Measured> measured;
    Size no_signal = 0, isotope_pattern = 0;

    // the unidentified multiplets, by the unique ids of their channels
    std::unordered_map<UInt64, const ConsensusFeature*> unidentified;
    for (const ConsensusFeature& cf : multiplets)
    {
      if (!cf.getPeptideIdentifications().empty()) { continue; }
      for (const FeatureHandle& handle : cf.getFeatures())
      {
        if (handle.hasValidUniqueId()) { unidentified.emplace(handle.getUniqueId(), &cf); }
      }
    }
    // Such a multiplet can be the candidate's peptide, detected without its identification mapping onto it, if it has its charge
    // and CV and its channels lie at the candidate's (detection tolerance; a complete one has one at the identified channel then).
    auto may_be = [&](const ConsensusFeature& cf, const TargetedCandidate_& candidate)
    {
      if (cf.getCharge() != candidate.charge || compensationVoltage_(cf) != candidate.cv) { return false; }
      return std::all_of(cf.getFeatures().begin(), cf.getFeatures().end(), [&](const FeatureHandle& handle)
                         {
                           return !handle.hasValidUniqueId()
                                  || std::any_of(candidate.mzs.begin(), candidate.mzs.end(), [&](double mz)
                                                 { return std::abs(handle.getMZ() - mz) <= detectionMZTolerance_(mz); });
                         });
    };

    for (auto& [key, group] : groups)
    {
      const std::string& peptide = key.first;
      std::stable_sort(group.begin(), group.end(), [](const TargetedCandidate_& a, const TargetedCandidate_& b) { return a.rt < b.rt; });

      // one feature per elution peak, with every identification in it (within the RT tolerance of its boundaries)
      vector<std::pair<TargetedPeak_, vector<Size>>> found; // the peak and its identifications (positions in the group)
      auto peak_at = [&](double rt)
      {
        return std::find_if(found.begin(), found.end(), [&](const auto& f) { return rt >= f.first.rt_lo - rt_tolerance && rt <= f.first.rt_hi + rt_tolerance; });
      };
      vector<TargetedOutcome_> outcome(group.size(), TargetedOutcome_::MEASURED);
      for (Size g = 0; g < group.size(); ++g)
      {
        if (auto it = peak_at(group[g].rt); it != found.end())
        {
          it->second.push_back(g);
          continue;
        }
        // Left out are the peaks of another identified peptide's pattern (OTHER_PEPTIDE, not isForeign(): see above), and those of
        // an unidentified multiplet that cannot be this peptide.
        const TargetedCandidate_& candidate = group[g];
        auto occupied = [&](const MSSpectrum& black, Size index)
        {
          if (peptides == nullptr) { return true; } // owners unknown: every reserved peak, as in the re-quantification
          const UInt64 owner = MultiplexResolverAlgorithm::blacklistOwner(black, index);
          if (MultiplexResolverAlgorithm::classifyBlacklistOwner(owner, peptide, *peptides) == MultiplexResolverAlgorithm::OwnerRelation::OTHER_PEPTIDE)
          {
            return true;
          }
          const auto it = unidentified.find(owner);
          return it != unidentified.end() && !may_be(*it->second, candidate);
        };
        TargetedPeak_ peak;
        outcome[g] = measureTargeted_(candidate, unassigned[candidate.index].getHits()[0], centroided, blacklist, occupied, elution, integrator, peak);
        if (outcome[g] != TargetedOutcome_::MEASURED) { continue; }
        auto it = std::find_if(found.begin(), found.end(), [&](const auto& f) { return peak.rt_lo <= f.first.rt_hi && peak.rt_hi >= f.first.rt_lo; });
        if (it != found.end()) { it->second.push_back(g); }
        else { found.emplace_back(std::move(peak), vector<Size>{g}); }
      }
      // an identification whose own channel gave no peak still belongs to a peak measured from another one
      for (Size g = 0; g < group.size(); ++g)
      {
        if (outcome[g] == TargetedOutcome_::MEASURED) { continue; }
        if (auto it = peak_at(group[g].rt); it != found.end())
        {
          it->second.push_back(g);
          outcome[g] = TargetedOutcome_::MEASURED;
        }
        else if (outcome[g] == TargetedOutcome_::NO_SIGNAL) { ++no_signal; }
        else { ++isotope_pattern; }
      }
      for (auto& [peak, members] : found) { measured.push_back({&group, std::move(peak), std::move(members)}); }
    }

    // One signal, one peptide: of the peaks of different peptides that read the same peaks (the identified channel of one at a
    // channel the other reads, same charge and CV, overlapping boundaries), the one with the best identification is kept.
    auto seed_of = [&](Size i) -> const TargetedCandidate_& { return (*measured[i].group)[measured[i].members.front()]; };
    vector<double> badness(measured.size(), std::numeric_limits<double>::infinity()); // the best score, lower = better
    for (Size i = 0; i < measured.size(); ++i)
    {
      for (Size g : measured[i].members)
      {
        const PeptideIdentification& id = unassigned[(*measured[i].group)[g].index];
        badness[i] = std::min(badness[i], id.isHigherScoreBetter() ? -id.getHits()[0].getScore() : id.getHits()[0].getScore());
      }
    }
    vector<Size> order(measured.size());
    std::iota(order.begin(), order.end(), Size(0));
    std::stable_sort(order.begin(), order.end(), [&](Size a, Size b) { return badness[a] < badness[b]; });
    auto same_signal = [&](Size a, Size b)
    {
      return measured[a].group != measured[b].group && seed_of(a).charge == seed_of(b).charge && seed_of(a).cv == seed_of(b).cv
             && measured[a].peak.rt_lo <= measured[b].peak.rt_hi && measured[a].peak.rt_hi >= measured[b].peak.rt_lo;
    };
    std::multimap<double, Size> kept_channels, kept_references; // m/z of the channels read and of the identified one, of the kept peaks
    auto shares = [&](const std::multimap<double, Size>& kept, double mz, Size i)
    {
      const double tolerance = detectionMZTolerance_(mz);
      for (auto it = kept.lower_bound(mz - tolerance); it != kept.end() && it->first <= mz + tolerance; ++it)
      {
        if (same_signal(i, it->second)) { return true; }
      }
      return false;
    };
    vector<bool> keep(measured.size(), false);
    Size other_identification = 0;
    for (Size i : order)
    {
      // a peak reads the peaks that end up in its feature: those of its measured channels (> 0) and of its identified channel (the
      // reference of the others, also if its own intensity is not quantifiable)
      const TargetedCandidate_& s = seed_of(i);
      const vector<double>& intensities = measured[i].peak.intensities;
      auto reads = [&](Size c) { return c == s.own || intensities[c] > 0.0; };
      bool conflict = shares(kept_channels, s.mzs[s.own], i);
      for (Size c = 0; !conflict && c < s.mzs.size(); ++c) { conflict = reads(c) && shares(kept_references, s.mzs[c], i); }
      if (conflict)
      {
        other_identification += measured[i].members.size(); // they stay unassigned
        continue;
      }
      keep[i] = true;
      for (Size c = 0; c < s.mzs.size(); ++c)
      {
        if (reads(c)) { kept_channels.emplace(s.mzs[c], i); }
      }
      kept_references.emplace(s.mzs[s.own], i);
    }

    Size used = 0, measured_channels = 0, absent_channels = 0, not_quantifiable = 0, identified_not_quantifiable = 0;
    vector<ConsensusFeature> features;
    vector<bool> moved(unassigned.size(), false);
    std::unordered_set<UInt64> read_owners; // features whose reserved peaks entered the targeted features
    for (Size i = 0; i < measured.size(); ++i)
    {
      if (!keep[i]) { continue; }
      const vector<TargetedCandidate_>& group = *measured[i].group;
      const TargetedPeak_& peak = measured[i].peak;
      vector<Size>& members = measured[i].members;
      const TargetedCandidate_& seed = group[members.front()]; // the identification the peak was measured for
      ConsensusFeature cf;
      for (Size c = 0; c < peak.intensities.size(); ++c)
      {
        const double intensity = peak.intensities[c];
        if (c != seed.own)
        {
          if (std::isnan(intensity)) { ++not_quantifiable; }
          else if (intensity > 0.0) { ++measured_channels; }
          else { ++absent_channels; }
        }
        else if (std::isnan(intensity)) { ++identified_not_quantifiable; }
        if (std::isnan(intensity)) { continue; }
        FeatureHandle handle;
        handle.setMapIndex(c);
        handle.setMZ(seed.mzs[c]);
        handle.setRT(peak.rt);
        handle.setIntensity(intensity);
        handle.setCharge(seed.charge);
        cf.insert(handle);
      }
      cf.setRT(peak.rt);
      cf.setMZ(seed.mzs.front());
      cf.setCharge(seed.charge);
      // the lightest measured channel is the anchor, as in a detected multiplet
      for (const FeatureHandle& handle : cf.getFeatures())
      {
        if (handle.getIntensity() > 0.0) { cf.setIntensity(handle.getIntensity()); break; }
      }
      if (seed.cv) { cf.setMetaValue(Constants::UserParam::FAIMS_CV, *seed.cv); }
      cf.setMetaValue("quantification_method", "targeted");
      std::sort(members.begin(), members.end());
      for (Size g : members)
      {
        PeptideIdentification& id = unassigned[group[g].index];
        id.setMetaValue("map_index", group[g].own);
        id.getHits()[0].setMetaValue("map_index", group[g].own);
        id.setMetaValue("quantification_method", "targeted");
        cf.getPeptideIdentifications().push_back(std::move(id));
        moved[group[g].index] = true;
        ++used;
      }
      cf.setUniqueId();
      features.push_back(std::move(cf));
      read_owners.insert(peak.owners.begin(), peak.owners.end());
    }

    // An unidentified multiplet kept for match between runs whose reserved peaks entered a targeted feature is that peptide (others'
    // peaks were left out): the feature replaces it, so that linking cannot give it the peptide's identification from another run.
    Size replaced = 0;
    if (!read_owners.empty())
    {
      const Size before = resolved.size();
      resolved.erase(std::remove_if(resolved.begin(), resolved.end(), [&](const ConsensusFeature& cf)
                                    {
                                      return cf.getPeptideIdentifications().empty()
                                             && std::any_of(cf.getFeatures().begin(), cf.getFeatures().end(), [&](const FeatureHandle& h)
                                                            { return h.hasValidUniqueId() && read_owners.contains(h.getUniqueId()); });
                                    }),
                     resolved.end());
      replaced = before - resolved.size();
    }

    if (!features.empty())
    {
      // a spectrum match belongs to one feature only: the moved identifications leave the unassigned ones
      PeptideIdentificationList kept;
      for (Size i = 0; i < unassigned.size(); ++i)
      {
        if (!moved[i]) { kept.push_back(std::move(unassigned[i])); }
      }
      unassigned = std::move(kept);
      for (ConsensusFeature& cf : features) { resolved.push_back(std::move(cf)); }
      resolved.sortByPosition();
    }
    OPENMS_LOG_INFO << "Targeted quantification: " << features.size() << " feature(s) for " << used << " identification(s) of peptides without a "
                    << "multiplet in this run (" << measured_channels << " other channel(s) measured against the identified one, " << absent_channels
                    << " absent, " << not_quantifiable << " not quantifiable); not quantified: " << already << " identification(s) of peptides quantified "
                    << "in this run, " << other_multiplet << " on a channel of another peptide's multiplet, " << other_identification
                    << " on the signal of a better identification of another peptide, " << no_signal << " without a clean elution peak of the "
                    << "identified channel, " << isotope_pattern << " whose isotope pattern does not fit, " << no_pattern << " whose labels fit no "
                    << "pattern, " << decoys << " decoy(s)." << endl;
    if (identified_not_quantifiable > 0)
    {
      OPENMS_LOG_INFO << "Of the targeted features, " << identified_not_quantifiable << " report their identified channel as not quantifiable: "
                      << "another channel's isotope makes up a substantial part of its peaks (the other channels are measured against them)." << endl;
    }
    if (replaced > 0)
    {
      OPENMS_LOG_INFO << replaced << " unidentified multiplet(s) kept for match between runs were replaced by the targeted feature(s) measured "
                      << "from their peaks." << endl;
    }
  }

  /// Let the peaks the detection reserved for a merged flank (blacklist with owners) belong to the channel it was merged into
  static void remapBlacklistOwners_(MSExperiment& blacklist, const std::unordered_map<UInt64, UInt64>& merged_into)
  {
    if (merged_into.empty() || !MultiplexResolverAlgorithm::hasBlacklistOwners(blacklist)) { return; }
    for (MSSpectrum& spectrum : blacklist)
    {
      std::vector<UInt64> owners(spectrum.size());
      bool changed = false;
      for (Size k = 0; k < spectrum.size(); ++k)
      {
        owners[k] = MultiplexResolverAlgorithm::blacklistOwner(spectrum, k);
        if (const auto it = merged_into.find(owners[k]); it != merged_into.end())
        {
          owners[k] = it->second;
          changed = true;
        }
      }
      if (changed) { MultiplexResolverAlgorithm::setBlacklistOwners(spectrum, owners); }
    }
  }

  /**
    @brief Quantify one run: multiplet detection, ID mapping, conflict resolution, multiplet completion.

    The result has one column per channel (map index = channel position), labelled with the channel's
    labels, and carries the run's identification run plus, as unassigned identifications, every
    identification that is not attached to a quantified multiplet.
  */
  ExitCodes quantifyRun_(const std::string& mz_file,
                         const std::string& id_file_abs_path,
                         const Size fraction_group,
                         const Size fraction,
                         const std::string& in_db,
                         set<std::string>& fixed_modifications,
                         set<std::string>& variable_modifications,
                         ConsensusMap& resolved)
  {
    OPENMS_LOG_INFO << "Processing " << File::basename(mz_file) << " (fraction group " << fraction_group << ", fraction " << fraction << ")" << endl;

    vector<ProteinIdentification> protein_ids;
    PeptideIdentificationList peptide_ids;
    MSExperiment ms1;
    ExitCodes e = loadAndCleanupIDFile_(id_file_abs_path, mz_file, in_db, fraction_group, fraction, protein_ids, peptide_ids,
                                        fixed_modifications, variable_modifications, ms1);
    if (e != EXECUTION_OK) { return e; }

    ConsensusMap multiplets;
    MSExperiment blacklist, centroided;
    FeatureMap channels;
    {
      FeatureMap features;
      e = detectMultiplets_(mz_file, std::move(ms1), multiplets, features, blacklist, centroided);
      if (e != EXECUTION_OK) { return e; }
      // partial multiplets that flank a larger one of the same peptide become part of it, with the peaks reserved for them;
      // every channel is described by where its multiplet elutes
      std::unordered_map<UInt64, UInt64> merged_into;
      channels = mergeFragments_(multiplets, features, merged_into);
      remapBlacklistOwners_(blacklist, merged_into);
    }

    // map the identifications onto the multiplets: an identification matches a multiplet if it lies in the elution
    // region of one of its channels (plus the tolerances), and the matched channel is recorded (needed by the resolver)
    {
      IDMapper mapper;
      Param idm_param = getParam_().copy("id_mapping:", true);
      writeDebug_("Parameters passed to IDMapper", idm_param, 3);
      mapper.setParameters(idm_param);
      MS1LabeledFAIMS::annotate(mapper, multiplets, channels, peptide_ids, protein_ids);
    }

    MultiplexResolverAlgorithm resolver;
    {
      Param res_param;
      res_param.insert("algorithm:", getParam_().copy("resolver:", true));
      res_param.setValue("algorithm:labels", labels_);
      res_param.setValue("algorithm:max_nr_labelled_aas", max_nr_labelled_aas_);
      // a missing channel that only the multiplet's own isotopes or dropped patterns occupy is left at 0 for the re-quantification
      res_param.setValue("algorithm:use_blacklist_owners", requantify_ ? "true" : "false");
      res_param.insert("labels:", getParam_().copy("label_mass_shifts:", true));
      writeDebug_("Parameters passed to MultiplexResolverAlgorithm", res_param, 3);
      resolver.setParameters(res_param);
    }

    if (knock_out_)
    {
      if (join_channels_)
      {
        // the channels of an identified peptide that were detected separately form one multiplet (with every spectrum match on
        // its closest multiplet only, so that a match the mapping attached to several does not join the wrong one)
        keepMultiAssignedOnClosest_(multiplets);
        joinChannelsOfIdentifiedPeptides_(multiplets, channels, resolver, centroided, blacklist);
      }
      // an identification of a channel the detection missed joins the multiplet detected in the peptide's other channels
      attachToUndetectedChannels_(multiplets, channels, resolver);
    }

    // a spectrum match that landed on several multiplets stays on the closest one only
    keepMultiAssignedOnClosest_(multiplets);

    // identifications whose labels fit the multiplet are preferred over better-scoring ones that do not
    keepConsistentIdentifications_(multiplets, resolver);

    // one identification (the best) per multiplet
    IDConflictResolverAlgorithm::resolve(multiplets);

    Size identified = 0;
    for (const auto& cf : multiplets) { identified += !cf.getPeptideIdentifications().empty(); }
    OPENMS_LOG_INFO << "Mapped identifications onto " << identified << " of " << multiplets.size() << " multiplet(s)." << endl;

    // remove quant/ID conflicts and complete the multiplets
    {
      ConsensusMap conflicts;
      resolver.resolve(multiplets, resolved, conflicts, blacklist);
      const Size n_resolved = resolved.size();

      Size conflicts_with_id = 0, unidentified = 0, kept_for_mbr = 0;
      for (auto& cf : conflicts)
      {
        if (cf.getPeptideIdentifications().empty())
        {
          ++unidentified;
          // With match between runs, a complete unidentified multiplet is kept: linking can hand it the
          // identification of a multiplet at the same position in another run. A complete multiplet was
          // detected by a complete pattern, whose channels are in the order of the columns; the channels
          // of a partial one are only known from an identification. (A channel added for an identification that
          // was then attached elsewhere has no unique id: it was not detected.)
          if (match_between_runs_
              && static_cast<Size>(std::count_if(cf.getFeatures().begin(), cf.getFeatures().end(),
                                                 [](const FeatureHandle& h) { return h.hasValidUniqueId(); })) == multiplicity_)
          {
            resolved.push_back(cf);
            ++kept_for_mbr;
          }
          continue;
        }
        // identifications of conflicting multiplets are still valid PSMs: keep them for inference, without quantification
        for (auto& id : cf.getPeptideIdentifications())
        {
          id.removeMetaValue("map_index");
          resolved.getUnassignedPeptideIdentifications().push_back(std::move(id));
          ++conflicts_with_id;
        }
      }
      if (kept_for_mbr > 0) { resolved.sortByPosition(); }
      OPENMS_LOG_INFO << "Resolved " << n_resolved << " multiplet(s); " << conflicts_with_id
                      << " identified multiplet(s) contradict their labels and are not quantified; " << unidentified
                      << " multiplet(s) are unidentified, " << kept_for_mbr << " of them kept for match between runs." << endl;
    }

    // the channels the detection did not find are measured where they are expected (without peaks reserved by other features,
    // as far as the blacklist tells their owners)
    if (requantify_ || targeted_quantification_)
    {
      const bool owners = MultiplexResolverAlgorithm::hasBlacklistOwners(blacklist);
      if (requantify_) { requantifyMissingChannels_(resolved, channels, centroided, blacklist, owners); }
      // identified peptides without a multiplet are measured where they are expected (while the handles still name their channels)
      if (targeted_quantification_)
      {
        const MultiplexResolverAlgorithm::OwnerPeptides peptides = owners ? resolver.ownerPeptides(multiplets) : MultiplexResolverAlgorithm::OwnerPeptides();
        quantifyTargeted_(resolved, multiplets, channels, resolver, centroided, blacklist, owners ? &peptides : nullptr);
      }
    }

    // The handles referred to the run's channel features by unique id; those are internal to this run and
    // not part of any output.
    for (auto& cf : resolved)
    {
      ConsensusFeature::HandleSetType handles;
      for (FeatureHandle handle : cf.getFeatures())
      {
        handle.setUniqueId(UniqueIdInterface::INVALID);
        handles.insert(handle);
      }
      cf.setFeatures(std::move(handles));
    }

    // A dummy channel whose region was blacklisted during feature detection is "not quantifiable" (NaN).
    // Report it as a missing channel value, like a run without a feature in label-free data: a NaN
    // would poison every abundance derived from it, and the QPX writers refuse non-finite values.
    {
      Size not_quantifiable = 0;
      for (auto& cf : resolved)
      {
        ConsensusFeature::HandleSetType kept;
        for (const FeatureHandle& handle : cf.getFeatures())
        {
          if (std::isnan(handle.getIntensity())) { ++not_quantifiable; continue; }
          kept.insert(handle);
        }
        if (kept.size() != cf.size()) { cf.setFeatures(std::move(kept)); }
      }
      if (not_quantifiable > 0)
      {
        OPENMS_LOG_INFO << not_quantifiable << " channel value(s) of completed multiplets are not quantifiable (another feature overlaps "
                        << "with the expected position, or the channel's peaks there do not form its isotope pattern) and are reported as missing." << endl;
      }
    }

    // The peptide identity is the unlabeled sequence: the label belongs to the channel, not to the
    // peptide, so the light and the heavy spectrum of one peptide must name one peptide for linking,
    // match between runs, inference and quantification. The resolver above has used the labels; from
    // here on every identification (quantified or not) is reduced to the peptide and carries its
    // label state as meta values.
    for (auto& cf : resolved)
    {
      for (auto& id : cf.getPeptideIdentifications())
      {
        for (auto& hit : id.getHits()) { stripLabels_(hit); }
      }
    }
    for (auto& id : resolved.getUnassignedPeptideIdentifications())
    {
      for (auto& hit : id.getHits()) { stripLabels_(hit); }
    }

    // the resolver records the channel of the identified feature of a completed multiplet on its hit;
    // make the identification's map index agree with it so that linking carries the right channel forward
    for (auto& cf : resolved)
    {
      for (auto& id : cf.getPeptideIdentifications())
      {
        if (!id.getHits().empty() && id.getHits()[0].metaValueExists("map_index"))
        {
          id.setMetaValue("map_index", id.getHits()[0].getMetaValue("map_index"));
        }
      }
    }

    // one column per channel of this run
    for (Size i = 0; i < multiplicity_; ++i)
    {
      ConsensusMap::ColumnHeader& header = resolved.getColumnHeaders()[i];
      header.filename = mz_file;
      header.label = channel_labels_[i];
      header.setMetaValue("channel_id", i);
      header.setMetaValue("channel_description", channel_descriptions_[i]); // the labels' PSI-MS names
      header.size = resolved.size();
    }
    resolved.setExperimentType("labeled_MS1");
    resolved.updateRanges();

    return EXECUTION_OK;
  }

  //-------------------------------------------------------------
  // per fraction: alignment and linking across runs
  //-------------------------------------------------------------

  /// Align the runs of one fraction by their identifications; returns the maximum RT difference left after alignment
  double align_(vector<ConsensusMap>& maps, vector<TransformationDescription>& transformations)
  {
    Param ma_param = getParam_().copy("alignment:", true);
    writeDebug_("Parameters passed to MapAlignmentAlgorithmIdentification", ma_param, 3);

    Param model_params = MapAlignerBase::getModelDefaults("b_spline");
    std::string model_type = model_params.getValue("type").toString();
    model_params = model_params.copy(model_type + ":", true);

    try
    {
      // the reference is determined from the data (the run that shares the most identifications with every other run)
      MapAlignmentAlgorithmIdentification aligner;
      aligner.setLogType(log_type_);
      aligner.setParameters(ma_param);
      aligner.align(maps, transformations, -1);
    }
    catch (Exception::MissingInformation& err)
    {
      if (getFlag_("force"))
      {
        OPENMS_LOG_ERROR << "Error: alignment failed. Details:\n" << err.what()
                         << "\nProcessing will continue using 'identity' transformations." << endl;
        model_type = "identity";
        transformations.clear();
        transformations.resize(maps.size());
      }
      else
      {
        throw;
      }
    }

    // fit the models (a NOP for 'identity'; too few points leave the transformation as identity as well)
    vector<TransformationDescription::TransformationStatistics> alignment_stats;
    for (TransformationDescription& t : transformations)
    {
      writeDebug_("Using " + StringUtils::toStr(t.getDataPoints().size()) + " points in fit.", 1);
      if (t.getDataPoints().size() > 10)
      {
        t.fitModel(model_type, model_params);
      }
      t.printSummary(getGlobalLogDebug());
      alignment_stats.emplace_back(t.getStatistics());
    }

    double max_alignment_diff = 0.0;
    for (const auto& s : alignment_stats)
    {
      if (s.percentiles_after.contains(100))
      {
        max_alignment_diff = std::max(max_alignment_diff, s.percentiles_after.at(100));
      }
    }
    OPENMS_LOG_INFO << "Max alignment difference (seconds): " << max_alignment_diff << endl;
    // very good alignments would otherwise leave no room for the linking; also cap the window
    max_alignment_diff = std::max(max_alignment_diff, 120.0);
    max_alignment_diff = std::min(max_alignment_diff, 600.0);
    return max_alignment_diff;
  }

  /// Link the runs of one fraction; the channels of every run are kept as sub-features (one column per run and channel)
  void link_(vector<ConsensusMap>& maps, double max_alignment_diff, ConsensusMap& consensus_fraction, const std::string& what = "consensus multiplet(s)")
  {
    Param fl_param = getParam_().copy("linking:", true);
    const double rt_typical = getParam_().getValue("algorithm:rt_typical");
    fl_param.setValue("distance_RT:max_difference", 2.0 * max_alignment_diff + 2.0 * rt_typical);
    writeDebug_("Parameters passed to FeatureGroupingAlgorithmQT", fl_param, 3);

    // the linker re-indexes the identifications by input map; remember the channel so that it can be restored below
    for (auto& map : maps)
    {
      map.applyFunctionOnPeptideIDs(
        [](PeptideIdentification& p)
        {
          if (p.metaValueExists("map_index")) { p.setMetaValue("old_map_index", p.getMetaValue("map_index")); }
        },
        true);
      map.updateRanges();
    }

    FeatureGroupingAlgorithmQT linker;
    linker.setParameters(fl_param);
    // All CVs of a physical run share its chromatographic RT transformation, but feature identity
    // and transfer between runs are restricted to the same compensation voltage.
    MS1LabeledFAIMS::group(linker, maps, consensus_fraction);

    consensus_fraction.applyMemberFunction(&UniqueIdInterface::setUniqueId);
    consensus_fraction.sortPeptideIdentificationsByMapIndex();
    OPENMS_LOG_INFO << "Linked " << maps.size() << " run(s) into " << consensus_fraction.size() << " " << what << "." << endl;
  }

  /**
    @brief Record on every multiplet of @p fraction the runs in which targeted quantification measured it

    The run names (file basenames) go to the meta value MS1LabeledRatioQuantifier::TARGETED_RUNS, which the abundances
    and ratios below read. After linking, a multiplet holds the channels of several runs, and only its identifications
    still tell which of them targeted quantification measured: those carry 'quantification_method' = 'targeted', and
    their 'map_index' names their column, i.e. their (run, channel). A targeted feature gets only such identifications,
    and a detected multiplet none of them. Targeted features are linked among themselves only (see alignAndLink_()), so a
    multiplet identified by nothing else holds targeted features only and counts as targeted in all of its runs.
  */
  static void markTargetedRuns_(ConsensusMap& fraction)
  {
    const auto& headers = fraction.getColumnHeaders();
    for (ConsensusFeature& cf : fraction)
    {
      std::set<std::string> runs;
      bool only_targeted = !cf.getPeptideIdentifications().empty();
      for (const PeptideIdentification& id : cf.getPeptideIdentifications())
      {
        if (id.getMetaValue("quantification_method", DataValue("")).toString() != "targeted")
        {
          only_targeted = false;
          continue;
        }
        if (!id.metaValueExists("map_index")) { continue; }
        const auto header = headers.find(static_cast<UInt64>(id.getMetaValue("map_index")));
        if (header != headers.end()) { runs.insert(File::basename(header->second.filename)); }
      }
      // If targeted features are all that identify the multiplet, it exists only because of targeted quantification, so all of
      // its runs count as targeted (a safeguard: they are linked among themselves only, so these are the runs of its IDs).
      if (only_targeted)
      {
        for (const FeatureHandle& handle : cf.getFeatures())
        {
          const auto header = headers.find(handle.getMapIndex());
          if (header != headers.end()) { runs.insert(File::basename(header->second.filename)); }
        }
      }
      if (!runs.empty()) { cf.setMetaValue(MS1LabeledRatioQuantifier::TARGETED_RUNS, StringList(runs.begin(), runs.end())); }
    }
  }

  /// Meta value (IntList) that requantifyMissingChannels_() puts on the identifications of a multiplet: the channels
  /// (0-based, the 'channel_id' of the run's columns) it re-quantified with a value; markRequantified_() moves it to the multiplet
  inline static const std::string requantified_channels_meta_ = "MS1Label:requantified_channel_ids";

  /**
    @brief Record on every multiplet of @p fraction the (run, channel) pairs that the re-quantification measured

    They go to the parallel meta values MS1LabeledRatioQuantifier::REQUANTIFIED_RUNS and REQUANTIFIED_LABELS, from which the
    ratios count how many of their evidences rest on a re-quantified channel (by run name, as the targeted runs: the columns
    change when the fractions are combined). After linking, only the identifications of a multiplet still tell its runs (their
    'map_index'); requantifyMissingChannels_() recorded the re-quantified channels on them. That record is removed from the
    identifications (also the unassigned ones).
  */
  static void markRequantified_(ConsensusMap& fraction)
  {
    const auto& headers = fraction.getColumnHeaders();
    for (ConsensusFeature& cf : fraction)
    {
      std::set<std::pair<std::string, Int>> requantified; // (run, label)
      for (PeptideIdentification& id : cf.getPeptideIdentifications())
      {
        if (!id.metaValueExists(requantified_channels_meta_)) { continue; }
        const IntList channels = id.getMetaValue(requantified_channels_meta_);
        id.removeMetaValue(requantified_channels_meta_);
        if (!id.metaValueExists("map_index")) { continue; }
        const auto header = headers.find(static_cast<UInt64>(id.getMetaValue("map_index")));
        if (header == headers.end()) { continue; }
        // the label of a channel is its channel_id + 1 (ConsensusMap::ColumnHeader::getLabelAsUInt())
        for (const Int channel : channels) { requantified.emplace(File::basename(header->second.filename), channel + 1); }
      }
      if (requantified.empty()) { continue; }
      StringList runs;
      IntList labels;
      for (const auto& [run, label] : requantified)
      {
        runs.push_back(run);
        labels.push_back(label);
      }
      cf.setMetaValue(MS1LabeledRatioQuantifier::REQUANTIFIED_RUNS, runs);
      cf.setMetaValue(MS1LabeledRatioQuantifier::REQUANTIFIED_LABELS, labels);
    }
    for (PeptideIdentification& id : fraction.getUnassignedPeptideIdentifications()) { id.removeMetaValue(requantified_channels_meta_); }
  }

  /**
    @brief Take the channels of targeted evidences that lack a value in some channel out of @p consensus

    A targeted feature whose other channel is not quantifiable keeps only its measured channel. Summed into the protein
    abundances, that channel would count without its partner and bias the abundance ratios (towards the identified
    channel, mostly at extreme ratios). Such evidences are therefore taken out before the abundances are read, and given
    back with restoreEvidences_() afterwards, so that the features keep them in the output.

    @return The handles taken out, by the index of their feature in @p consensus
  */
  std::vector<std::pair<Size, std::vector<FeatureHandle>>> withholdIncompleteTargeted_(ConsensusMap& consensus) const
  {
    std::vector<std::pair<Size, std::vector<FeatureHandle>>> withheld;
    Size evidences = 0;
    const auto& headers = consensus.getColumnHeaders();
    for (Size i = 0; i < consensus.size(); ++i)
    {
      ConsensusFeature& cf = consensus[i];
      if (!cf.metaValueExists(MS1LabeledRatioQuantifier::TARGETED_RUNS)) { continue; }
      const StringList targeted = cf.getMetaValue(MS1LabeledRatioQuantifier::TARGETED_RUNS);
      const std::set<std::string> targeted_runs(targeted.begin(), targeted.end());

      // the handles of every targeted run of this feature
      std::map<std::string, std::vector<FeatureHandle>> by_run;
      for (const FeatureHandle& handle : cf.getFeatures())
      {
        const auto header = headers.find(handle.getMapIndex());
        if (header == headers.end()) { continue; }
        const std::string run = File::basename(header->second.filename);
        if (targeted_runs.contains(run)) { by_run[run].push_back(handle); }
      }
      std::vector<FeatureHandle> taken;
      std::set<UInt64> taken_columns;
      for (auto& [run, handles] : by_run)
      {
        // a not-quantifiable channel has no handle (see quantifyRun_), an absent one has intensity 0
        if (handles.size() >= multiplicity_) { continue; }
        ++evidences;
        for (const FeatureHandle& handle : handles)
        {
          taken_columns.insert(handle.getMapIndex());
          taken.push_back(handle);
        }
      }
      if (taken.empty()) { continue; }
      ConsensusFeature::HandleSetType kept;
      for (const FeatureHandle& handle : cf.getFeatures())
      {
        if (!taken_columns.contains(handle.getMapIndex())) { kept.insert(handle); }
      }
      cf.setFeatures(std::move(kept));
      withheld.emplace_back(i, std::move(taken));
    }
    if (evidences > 0)
    {
      OPENMS_LOG_INFO << "Abundances: " << evidences << " evidence(s) of targeted quantification with a channel that is not quantifiable "
                      << "are left out of the peptide and protein abundances." << endl;
    }
    return withheld;
  }

  /// Give the handles that withholdIncompleteTargeted_() took out back to their features
  static void restoreEvidences_(ConsensusMap& consensus, std::vector<std::pair<Size, std::vector<FeatureHandle>>>& withheld)
  {
    for (auto& [index, handles] : withheld)
    {
      for (const FeatureHandle& handle : handles) { consensus[index].insert(handle); }
    }
    withheld.clear();
  }

  /**
    @brief Move the features of targeted quantification out of @p maps

    @return One map per run with the targeted features of that run, the same columns and protein runs and no unassigned identifications
  */
  static vector<ConsensusMap> splitOffTargeted_(vector<ConsensusMap>& maps)
  {
    vector<ConsensusMap> targeted(maps.size());
    for (Size i = 0; i < maps.size(); ++i)
    {
      ConsensusMap& map = maps[i];
      ConsensusMap& part = targeted[i];
      part = map;
      part.resize(0);
      part.getUnassignedPeptideIdentifications().clear();
      for (auto& [index, header] : part.getColumnHeaders()) { header.size = 0; }
      ConsensusMap::iterator kept = map.begin();
      for (ConsensusMap::iterator it = map.begin(); it != map.end(); ++it)
      {
        if (it->getMetaValue("quantification_method", DataValue("")).toString() == "targeted")
        {
          for (const FeatureHandle& handle : it->getFeatures())
          {
            ++part.getColumnHeaders().at(handle.getMapIndex()).size;
            --map.getColumnHeaders().at(handle.getMapIndex()).size;
          }
          part.push_back(std::move(*it));
        }
        else
        {
          if (kept != it) { *kept = std::move(*it); }
          ++kept;
        }
      }
      map.erase(kept, map.end());
      map.updateRanges();
      part.updateRanges();
    }
    return targeted;
  }

  /**
    @brief Align and link the runs of one fraction

    The features of targeted quantification take no part in the alignment and are linked among themselves only, after the
    detected multiplets have been linked. So the alignment, the linking of the detected multiplets and with match between
    runs the identifications that their unidentified multiplets take over are those without targeted quantification:
    a targeted feature, whose position rests on one identification, cannot widen the RT tolerance of the linking or draw an
    unidentified multiplet of another peptide into its linked feature.
  */
  void alignAndLink_(vector<ConsensusMap>& maps, ConsensusMap& consensus_fraction)
  {
    if (maps.size() > 1)
    {
      vector<ConsensusMap> targeted;
      if (targeted_quantification_) { targeted = splitOffTargeted_(maps); }
      const bool any_targeted = std::any_of(targeted.begin(), targeted.end(), [](const ConsensusMap& m) { return !m.empty(); });

      vector<TransformationDescription> transformations;
      const double max_alignment_diff = align_(maps, transformations);
      const auto transform = [](ConsensusMap& map, const TransformationDescription& transformation)
      {
        try
        {
          MapAlignmentTransformer::transformRetentionTimes(map, transformation);
        }
        catch (Exception::IllegalArgument& e)
        {
          OPENMS_LOG_WARN << e.what() << endl;
        }
      };
      for (Size i = 0; i < maps.size(); ++i)
      {
        transform(maps[i], transformations[i]);
        if (any_targeted) { transform(targeted[i], transformations[i]); }
      }
      link_(maps, max_alignment_diff, consensus_fraction);
      if (any_targeted)
      {
        ConsensusMap linked_targeted;
        link_(targeted, max_alignment_diff, linked_targeted, "consensus feature(s) of targeted quantification");
        MS1LabeledFAIMS::append(consensus_fraction, std::move(linked_targeted));
        consensus_fraction.updateRanges();
      }
    }
    else
    {
      // a fraction measured in a single run: nothing to align or link
      consensus_fraction = std::move(maps.front());
    }
  }

  //-------------------------------------------------------------
  // protein inference and FDR (see ProteomicsLFQ)
  //-------------------------------------------------------------

  ExitCodes inferProteinGroups_(ConsensusMap& consensus, const set<std::string>& fixed_modifications)
  {
    const bool groups = getStringOption_("protein_quantification") != "strictly_unique_peptides";
    const bool bayesian = getStringOption_("protein_inference") == "bayesian";
    const bool greedy_group_resolution = getStringOption_("protein_quantification") == "shared_peptides";

    // study-wide inference operates on a single merged ID run
    ConsensusMapMergerAlgorithm cmerge;
    cmerge.mergeAllIDRuns(consensus);

    if (!bayesian)
    {
      BasicProteinInferenceAlgorithm bpia;
      Param bpiaparams = bpia.getParameters();
      bpiaparams.setValue("annotate_indistinguishable_groups", groups ? "true" : "false");
      bpiaparams.setValue("greedy_group_resolution", greedy_group_resolution ? "true" : "false");
      bpia.setParameters(bpiaparams);
      bpia.run(consensus, consensus.getProteinIdentifications()[0], true);
    }
    else
    {
      BayesianProteinInferenceAlgorithm bayes;
      Param bayesparams = bayes.getParameters();
      bayesparams.setValue("keep_best_PSM_only", "false");
      bayes.setParameters(bayesparams);
      bayes.inferPosteriorProbabilities(consensus, greedy_group_resolution);
      if (!groups)
      {
        consensus.getProteinIdentifications()[0].getIndistinguishableProteins().clear();
      }
    }

    const double max_fdr = getDoubleOption_("proteinFDR");
    const bool picked = getStringOption_("picked_proteinFDR") == "true";
    const double max_psm_fdr = getDoubleOption_("psmFDR");
    FalseDiscoveryRate fdr;

    auto& overall_proteins = consensus.getProteinIdentifications()[0];
    if (!picked) { fdr.applyBasic(overall_proteins); }
    else { fdr.applyPickedProteinFDR(overall_proteins, picked_decoy_string_, picked_decoy_prefix_); }

    const bool pepFDR = getStringOption_("FDR_type") == "PSM+peptide";
    if (pepFDR) { fdr.applyBasicPeptideLevel(consensus, true); }
    else { fdr.applyBasic(consensus, true); }

    // FDR filtering removed all decoy proteins -> update references and remove all unreferenced (decoy) PSMs
    IDFilter::removeDanglingProteinReferences(consensus, true);
    IDFilter::removeUnreferencedProteins(consensus, true);
    IDFilter::updateProteinGroups(overall_proteins.getIndistinguishableProteins(), overall_proteins.getHits());
    IDFilter::updateProteinGroups(overall_proteins.getProteinGroups(), overall_proteins.getHits());

    if (max_psm_fdr < 1.0)
    {
      for (auto& f : consensus)
      {
        IDFilter::filterHitsByScore(f.getPeptideIdentifications(), max_psm_fdr);
      }
      IDFilter::filterHitsByScore(consensus.getUnassignedPeptideIdentifications(), max_psm_fdr);
    }
    if (max_fdr < 1.0)
    {
      IDFilter::filterHitsByScore(overall_proteins, max_fdr);
    }
    IDFilter::removeDanglingProteinReferences(consensus, true);
    if (max_psm_fdr < 1.0) { IDFilter::removeUnreferencedProteins(consensus, true); }
    IDFilter::updateProteinGroups(overall_proteins.getIndistinguishableProteins(), overall_proteins.getHits());
    IDFilter::updateProteinGroups(overall_proteins.getProteinGroups(), overall_proteins.getHits());

    if (overall_proteins.getHits().empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "No proteins left after FDR filtering. Please check the log and adjust your settings.");
    }

    // strictly unique peptides: filter for the theoretical uniqueness annotated during peptide indexing
    if (!greedy_group_resolution && !groups)
    {
      for (auto& f : consensus)
      {
        IDFilter::keepUniquePeptidesPerProtein(f.getPeptideIdentifications());
      }
      IDFilter::keepUniquePeptidesPerProtein(consensus.getUnassignedPeptideIdentifications());
      IDFilter::removeUnreferencedProteins(consensus, true);

      if (overall_proteins.getHits().empty())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "No protein is supported by a strictly unique peptide. Note that "
                                            "'protein_quantification' = 'strictly_unique_peptides' requires the "
                                            "theoretical uniqueness annotated during peptide indexing (see '-fasta').");
      }
      // quantification and every exporter report abundances on protein groups
      overall_proteins.fillIndistinguishableGroupsWithSingletons();
    }

    // coverage needs the protein sequences, which only peptide indexing ('-fasta') or an indexed input provides
    const bool sequences_known = std::all_of(overall_proteins.getHits().begin(), overall_proteins.getHits().end(),
                                             [](const ProteinHit& h) { return !h.getSequence().empty(); });
    if (sequences_known)
    {
      overall_proteins.computeCoverage(consensus, true);
    }
    else
    {
      OPENMS_LOG_INFO << "Protein sequences are not available (no '-fasta'), skipping the coverage computation." << endl;
    }

    // determine observed modifications (exclude fixed mods)
    overall_proteins.computeModifications(consensus, StringList(fixed_modifications.begin(), fixed_modifications.end()), true);

    return EXECUTION_OK;
  }

  //-------------------------------------------------------------
  // label handling
  //-------------------------------------------------------------

  /// Channel label as FeatureFinderMultiplex writes it into the column header, e.g. "no_label" or "Lys8Arg10"
  static std::string channelLabel_(const vector<std::string>& sample_labels)
  {
    std::string label;
    for (const std::string& l : sample_labels) { label += l; }
    return label;
  }

  /// Derive the channels and the modifications the search has to contain from '-labels'
  void initLabels_()
  {
    labels_ = getStringOption_("labels");
    max_nr_labelled_aas_ = getIntOption_("max_nr_labelled_aas");

    std::map<std::string, double> label_mass_shift;
    Param shifts = getParam_().copy("label_mass_shifts:", true);
    for (Param::ParamIterator it = shifts.begin(); it != shifts.end(); ++it)
    {
      label_mass_shift[it->name] = static_cast<double>(it->value);
    }

    // throws Exception::InvalidParameter for unknown labels
    MultiplexDeltaMassesGenerator generator(labels_, max_nr_labelled_aas_, label_mass_shift);
    const vector<vector<std::string>> samples = generator.getSamplesLabelsList();

    multiplicity_ = samples.size();
    if (multiplicity_ < 2)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "'-labels " + labels_ + "' describes a single, unlabelled channel. This tool quantifies MS1-labeled experiments with at "
        "least two channels; use ProteomicsLFQ for label-free data.");
    }

    // The generator may insert an implicit light SILAC channel and omit channels that cannot
    // occur in a pattern. Every detected handle must still have the channel promised by the design.
    if (generator.getDeltaMassesList().empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "'-labels " + labels_ + "' does not generate any mass patterns. Check the label mass shifts.");
    }
    for (const auto& pattern : generator.getDeltaMassesList())
    {
      const auto& masses = pattern.getDeltaMasses();
      bool consistent = masses.size() == samples.size();
      for (Size i = 0; consistent && i < masses.size(); ++i)
      {
        for (const auto& label : masses[i].label_set)
        {
          if (std::find(samples[i].begin(), samples[i].end(), label) == samples[i].end()) { consistent = false; }
        }
        if (i > 0 && masses[i].delta_mass <= masses[i - 1].delta_mass) { consistent = false; }
      }
      if (! consistent)
      {
        throw Exception::InvalidParameter(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "The channel order in '-labels " + labels_
            + "' does not match the generated mass patterns. "
              "Specify every channel in increasing mass order, with an explicit unlabelled first channel [] for SILAC "
              "(e.g. [][Lys4,Arg6][Lys8,Arg10]). Each pattern must contain every declared channel.");
      }
    }

    channel_labels_.clear();
    channel_label_sets_.clear();
    channel_descriptions_.clear();
    required_label_mods_.clear();
    label_long_to_shorts_.clear();
    set<std::string> seen;
    for (const auto& sample : samples)
    {
      channel_labels_.push_back(channelLabel_(sample));
      channel_label_sets_.emplace_back(sample.begin(), sample.end());
      vector<std::string> description;
      for (const std::string& short_name : sample)
      {
        if (short_name == "no_label")
        {
          description.push_back("unlabeled");
          continue;
        }
        const std::string long_name = generator.getLabelLong(short_name);
        // a plain mass shift has no modification name to check for or to strip
        if (long_name.empty() || long_name == short_name)
        {
          description.push_back(short_name);
          continue;
        }
        description.push_back(long_name);
        if (seen.insert(short_name).second)
        {
          required_label_mods_.emplace_back(short_name, long_name);
          label_long_to_shorts_[long_name].push_back(short_name);
        }
      }
      channel_descriptions_.push_back(ListUtils::concatenate(description, ", "));
    }
    OPENMS_LOG_INFO << "Labels: " << multiplicity_ << " channels (" << ListUtils::concatenate(channel_labels_, ", ") << ")" << endl;
  }

  /// The short name of a label modification on @p residue ('\0' for a terminus), or "" if @p mod is not one of '-labels'
  std::string labelShortName_(const ResidueModification& mod, char residue) const
  {
    const auto it = label_long_to_shorts_.find(mod.getId());
    if (it == label_long_to_shorts_.end()) { return ""; }
    // Arg6 and Lys6 share one modification (Label:13C(6)); the residue tells them apart
    for (const std::string& short_name : it->second)
    {
      if ((residue == 'K' && StringUtils::hasPrefix(short_name, "Lys")) || (residue == 'R' && StringUtils::hasPrefix(short_name, "Arg")))
      {
        return short_name;
      }
    }
    return it->second.front();
  }

  /// The 1-based channel whose labels cover @p removed (the labels a peptidoform carried); 0 if none or several do
  int channelOfLabels_(const vector<std::string>& removed) const
  {
    set<std::string> labels(removed.begin(), removed.end());
    if (labels.empty()) { labels.insert("no_label"); }
    int channel = 0;
    for (Size i = 0; i < channel_label_sets_.size(); ++i)
    {
      if (std::includes(channel_label_sets_[i].begin(), channel_label_sets_[i].end(), labels.begin(), labels.end()))
      {
        if (channel != 0) { return 0; } // ambiguous
        channel = static_cast<int>(i + 1);
      }
    }
    return channel;
  }

  /// The sequence of @p hit without its labels (other modifications kept), see stripLabels_
  AASequence unlabelled_(const PeptideHit& hit) const
  {
    PeptideHit copy = hit;
    stripLabels_(copy);
    return copy.getSequence();
  }

  /**
    @brief Reduce the hit's sequence to the peptide identity and record its label state.

    The label modifications of '-labels' are removed from the sequence (the channel, not the
    peptide, carries the label), and the hit is annotated with 'labeled_sequence' (as searched),
    'removed_labels' (short names, or 'none') and 'channel' (1-based; the 'Label' of the
    experimental design; 0 if unknown). These are exported as opt_ columns by mzTab and as
    cv_params by the QPX feature and psm views.
  */
  void stripLabels_(PeptideHit& hit) const
  {
    const AASequence original = hit.getSequence();
    AASequence stripped = original;
    vector<std::string> removed;

    if (stripped.hasNTerminalModification())
    {
      const std::string short_name = labelShortName_(*stripped.getNTerminalModification(), '\0');
      if (!short_name.empty())
      {
        removed.push_back(short_name);
        stripped.setNTerminalModification("");
      }
    }
    for (Size i = 0; i < stripped.size(); ++i)
    {
      const Residue& residue = stripped[i];
      if (!residue.isModified()) { continue; }
      const std::string short_name = labelShortName_(*residue.getModification(), residue.getOneLetterCode()[0]);
      if (!short_name.empty())
      {
        removed.push_back(short_name);
        stripped.setModification(i, "");
      }
    }
    if (stripped.hasCTerminalModification())
    {
      const std::string short_name = labelShortName_(*stripped.getCTerminalModification(), '\0');
      if (!short_name.empty())
      {
        removed.push_back(short_name);
        stripped.setCTerminalModification("");
      }
    }

    hit.setSequence(stripped);
    hit.setMetaValue(MS1LabelState::LABELED_SEQUENCE, original.toString());
    hit.setMetaValue(MS1LabelState::REMOVED_LABELS, removed.empty() ? std::string("none") : ListUtils::concatenate(removed, ","));
    hit.setMetaValue(MS1LabelState::CHANNEL, channelOfLabels_(removed));
  }

  /// QPX names channels with a fixed vocabulary; refuse '-out_qpx' up front for labels outside of it
  void checkQPXLabels_() const
  {
    ConsensusMap probe;
    probe.setExperimentType("labeled_MS1");
    for (Size i = 0; i < multiplicity_; ++i)
    {
      ConsensusMap::ColumnHeader& header = probe.getColumnHeaders()[i];
      header.filename = "probe.mzML";
      header.label = channel_labels_[i];
      header.setMetaValue("channel_id", i);
    }
    const auto qpx_labels = ArrowIOHelpers::qpxIntensityLabels(probe);
    for (Size i = 0; i < multiplicity_; ++i)
    {
      if (qpx_labels.at(i).empty())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "'-out_qpx' cannot be written for '-labels " + labels_ + "': channel " + StringUtils::toStr(i + 1) + " ('" + channel_labels_[i]
          + "') has no label in the QPX intensity-label vocabulary (SILAC light/medium/heavy for two- and three-plex SILAC, "
          "DIMETHYL0/2/4/6/8). Use '-out' or '-out_cxml' instead.");
      }
    }
  }

  //-------------------------------------------------------------
  // main
  //-------------------------------------------------------------

  ExitCodes main_(int, const char**) override
  {
    const std::string out = getStringOption_("out");
    const std::string out_cxml = getStringOption_("out_cxml");
    const std::string out_qpx = getOutputDirOption("out_qpx");
    if (out.empty() && out_cxml.empty() && out_qpx.empty())
    {
      throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "out/out_cxml/out_qpx");
    }

    StringList in = getStringList_("in");
    StringList in_ids = getStringList_("ids");
    if (in.size() != in_ids.size())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Number of spectra files (" + StringUtils::toStr(in.size()) + ") must match number of ID files (" + StringUtils::toStr(in_ids.size()) + ").");
    }
    // Design lookup and downstream exporters identify physical runs by basename. Reject aliases
    // before any design filtering or path reconciliation can silently collapse distinct inputs.
    set<std::string> in_basenames;
    for (const auto& file : in)
    {
      if (! in_basenames.insert(File::basename(file)).second)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Input spectra files must have unique basenames; duplicate basename '" + File::basename(file)
                                            + "'. Rename the files and update the experimental design before processing.");
      }
    }
    const std::string design_file = getStringOption_("design");
    const std::string in_db = getStringOption_("fasta");

    initLabels_();
    if (!out_qpx.empty()) { checkQPXLabels_(); }

    if (const int reference_channel = getParam_().getValue("ratios:reference_channel"); reference_channel > static_cast<int>(multiplicity_))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "'-ratios:reference_channel " + StringUtils::toStr(reference_channel) + "' names a channel that '-labels " + labels_ + "' does not have ("
        + StringUtils::toStr(multiplicity_) + " channels).");
    }

    match_between_runs_ = getStringOption_("match_between_runs") == "true";
    requantify_ = getStringOption_("requantify") == "true";
    join_channels_ = getStringOption_("join_channels") == "true";
    requantify_min_correlation_ = getDoubleOption_("requantify_min_correlation");
    requantify_correlation_ratio_ = getDoubleOption_("requantify_correlation_ratio");
    targeted_quantification_ = getStringOption_("targeted_quantification") == "true";
    knock_out_ = getParam_().getValue("algorithm:knock_out").toString() == "true";
    if (knock_out_ && multiplicity_ > 4)
    {
      OPENMS_LOG_WARN << "Warning: the detection of partial multiplets ('algorithm:knock_out') supports up to four channels; '-labels "
                      << labels_ << "' has " << multiplicity_ << ". Only multiplets in which every channel is detected are quantified." << endl;
      knock_out_ = false;
    }

    if (getStringOption_("protein_quantification") == "strictly_unique_peptides" && in_db.empty())
    {
      OPENMS_LOG_WARN << "Warning: '-protein_quantification strictly_unique_peptides' filters peptides by their theoretical uniqueness, "
                         "which is annotated during peptide indexing. Without '-fasta' this relies on the input identifications already "
                         "being indexed (e.g., by PeptideIndexer); otherwise no peptide will pass the filter." << endl;
    }

    //-------------------------------------------------------------
    // experimental design: read or generate default
    //-------------------------------------------------------------
    if (!design_file.empty())
    {
      design_ = ExperimentalDesignFile::load(design_file, false);
    }
    else
    {
      OPENMS_LOG_INFO << "No experimental design file provided.\n"
                      << "Assuming an unfractionated experiment: one fraction group per file, one sample per (file, channel).\n" << endl;
      TextFile design_table;
      design_table.addLine("Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample");
      Size sample = 1;
      for (Size i = 0; i < in.size(); ++i)
      {
        for (Size label = 1; label <= multiplicity_; ++label)
        {
          design_table.addLine(StringUtils::toStr(i + 1) + "\t1\t" + in[i] + "\t" + StringUtils::toStr(label) + "\t" + StringUtils::toStr(sample++));
        }
      }
      design_ = ExperimentalDesignFile::load(design_table, false, "--no design file--");
    }

    // A basename must also identify a single path in the selected design. Otherwise two design
    // runs could both be reconciled to the same input file even when the inputs themselves are unique.
    std::map<std::string, std::string> design_paths;
    for (const auto& entry : design_.getMSFileSection())
    {
      const auto basename = File::basename(entry.path);
      if (! in_basenames.contains(basename)) { continue; }
      const auto [it, inserted] = design_paths.emplace(basename, entry.path);
      if (! inserted && it->second != entry.path)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "The experimental design contains distinct paths with the same basename '" + basename
                                            + "'. Use unique basenames for physical runs.");
      }
    }

    // the design's channels must be the ones of '-labels'
    if (design_.getNumberOfLabels() != multiplicity_)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "The experimental design uses " + StringUtils::toStr(design_.getNumberOfLabels()) + " label(s), but '-labels " + labels_ + "' describes "
        + StringUtils::toStr(multiplicity_) + " channels. The 'Label' column has to number the channels of '-labels' from 1.");
    }

    // input files must be covered by the design, with every channel
    {
      const auto& pl2fg = design_.getPathLabelToFractionGroupMapping(true);
      for (const std::string& bn : in_basenames)
      {
        for (Size label = 1; label <= multiplicity_; ++label)
        {
          if (!pl2fg.contains({bn, static_cast<unsigned>(label)}))
          {
            throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "The experimental design has no row for spectra file '" + bn + "' and label " + StringUtils::toStr(label)
              + ". Every input file needs one row per channel of '-labels'.");
          }
        }
      }
    }
    const Size nr_filtered = design_.filterByBasenames(in_basenames);
    if (nr_filtered > 0)
    {
      OPENMS_LOG_WARN << "WARNING: " << nr_filtered
                      << " files from experimental design were not passed as spectra inputs. Continuing with subset if the fractions still match." << endl;
    }
    if (!design_.sameNrOfMSFilesPerFraction())
    {
      OPENMS_LOG_WARN << "WARNING: Different number of fractions for different samples provided. Support maybe limited." << endl;
    }

    std::map<unsigned int, std::vector<std::string>> frac2ms = design_.getFractionToMSFilesMapping();
    for (auto& fraction_ms_files : frac2ms)
    {
      // a multiplexed file has one design row per channel; it is still one run
      std::vector<std::string>& files = fraction_ms_files.second;
      std::vector<std::string> unique_files;
      for (const std::string& f : files)
      {
        if (std::find(unique_files.begin(), unique_files.end(), f) == unique_files.end()) { unique_files.push_back(f); }
      }
      files.swap(unique_files);

      // the design may name the files with other paths; use the ones that exist
      for (auto& s : fraction_ms_files.second)
      {
        if (auto it = std::find_if(in.begin(), in.end(), [&s](const std::string& in_filename) { return File::basename(in_filename) == File::basename(s); });
            it != in.end() && s != *it)
        {
          OPENMS_LOG_INFO << "Path of spectra files differ between experimental design (1) and input (2). Using the path of the input file as "
                          << "we know this file exists on the file system: '" << *it << "' vs. '" << s << endl;
          s = *it;
        }
      }
    }

    // spectra file -> ID file (same order)
    const std::map<std::string, std::string> mzfile2idfile = DDAWorkflowCommons::mapMzML2Ids(in, in_ids);

    //-------------------------------------------------------------
    // quantification per fraction
    //-------------------------------------------------------------
    set<std::string> fixed_modifications, variable_modifications;
    ConsensusMap consensus;
    const auto& path_label_to_fractiongroup = design_.getPathLabelToFractionGroupMapping(true);

    for (const auto& [fraction, ms_files] : frac2ms)
    {
      OPENMS_LOG_INFO << "Fraction " << fraction << ": " << ms_files.size() << " run(s)" << endl;

      vector<ConsensusMap> maps;
      maps.reserve(ms_files.size());
      for (const std::string& mz_file : ms_files)
      {
        const Size fraction_group = path_label_to_fractiongroup.at({File::basename(mz_file), 1});
        const std::string id_file_abs_path = File::absolutePath(mzfile2idfile.at(File::absolutePath(mz_file)));

        ConsensusMap resolved;
        ExitCodes e = quantifyRun_(mz_file, id_file_abs_path, fraction_group, fraction, in_db, fixed_modifications, variable_modifications, resolved);
        if (e != EXECUTION_OK) { return e; }
        maps.push_back(std::move(resolved));
      }

      ConsensusMap consensus_fraction;
      alignAndLink_(maps, consensus_fraction);
      if (targeted_quantification_) { markTargetedRuns_(consensus_fraction); }
      if (requantify_) { markRequantified_(consensus_fraction); }

      // fractions are combined column-wise: one column per (run, channel), rows stay per fraction
      consensus.appendColumns(consensus_fraction);
    }

    consensus.setExperimentType("labeled_MS1");
    {
      StringList ms_runs;
      for (const auto& [map_index, header] : consensus.getColumnHeaders()) { ms_runs.push_back(header.filename); }
      consensus.setPrimaryMSRunPath(ms_runs);
    }
    consensus.sortByPosition();
    consensus.sortPeptideIdentificationsByMapIndex();

    // the fraction structure lives in the design; put it on the columns for the exporters
    if (const Size unannotated = design_.annotateColumnHeaders(consensus); unannotated > 0)
    {
      OPENMS_LOG_WARN << unannotated << " of " << consensus.getColumnHeaders().size()
                      << " (file, channel) column(s) have no matching row in the experimental design." << endl;
    }

    if (consensus.empty())
    {
      OPENMS_LOG_FATAL_ERROR << "No identified multiplet was quantified in any run. Check that the labels were searched as modifications "
                                "and that the identifications and spectra belong together." << endl;
      return UNEXPECTED_RESULT;
    }

    //-------------------------------------------------------------
    // protein inference and quantification
    //-------------------------------------------------------------
    ExitCodes e = inferProteinGroups_(consensus, fixed_modifications);
    if (e != EXECUTION_OK) { return e; }

    IDFilter::removeUnreferencedProteins(consensus, true);
    IDConflictResolverAlgorithm::resolve(consensus);

    // The label state of the identification a multiplet is quantified under, at feature level as well:
    // consensusXML and the mzTab peptide section report feature meta values, the hit's own stay with
    // the PSM. After linking this identification may come from another run (match between runs).
    for (auto& cf : consensus)
    {
      if (cf.getPeptideIdentifications().empty() || cf.getPeptideIdentifications()[0].getHits().empty()) { continue; }
      const PeptideHit& hit = cf.getPeptideIdentifications()[0].getHits()[0];
      for (const std::string& key : MS1LabelState::keys())
      {
        if (hit.metaValueExists(key)) { cf.setMetaValue(key, hit.getMetaValue(key)); }
      }
    }

    Param pq_param = getParam_().copy("ProteinQuantification:", true);
    writeDebug_("Parameters passed to PeptideAndProteinQuant algorithm", pq_param, 3);
    PeptideAndProteinQuant quantifier;
    quantifier.setParameters(pq_param);
    {
      auto withheld = withholdIncompleteTargeted_(consensus);
      quantifier.readQuantData(consensus, design_);
      restoreEvidences_(consensus, withheld);
    }
    quantifier.quantifyPeptides();

    ProteinIdentification& inferred_proteins = consensus.getProteinIdentifications()[0];
    if (inferred_proteins.getIndistinguishableProteins().empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No information on indistinguishable protein groups found.");
    }
    quantifier.quantifyProteins(inferred_proteins);
    const auto& protein_quants = quantifier.getProteinResults();
    if (protein_quants.empty()) { OPENMS_LOG_WARN << "Warning: No proteins were quantified." << endl; }
    // keep protein groups that have not been quantified, as the sibling workflows do
    quantifier.annotateQuantificationsToProteins(protein_quants, inferred_proteins, false);

    // The quantity of a labeled experiment is the channel ratio, aggregated as a median of ratios
    // rather than as a ratio of aggregated intensities (see MS1LabeledRatioQuantifier). Runs after
    // the abundances above, so that both sit next to each other on the same protein groups.
    {
      MS1LabeledRatioQuantifier ratio_quantifier;
      Param ratio_param = getParam_().copy("ratios:", true);
      writeDebug_("Parameters passed to MS1LabeledRatioQuantifier", ratio_param, 3);
      ratio_quantifier.setParameters(ratio_param);
      ratio_quantifier.run(consensus, design_, inferred_proteins);
    }
    // the targeted runs and re-quantified columns were bookkeeping for the abundances and ratios; the identifications keep 'quantification_method'
    for (auto& cf : consensus)
    {
      cf.removeMetaValue(MS1LabeledRatioQuantifier::TARGETED_RUNS);
      cf.removeMetaValue(MS1LabeledRatioQuantifier::REQUANTIFIED_RUNS);
      cf.removeMetaValue(MS1LabeledRatioQuantifier::REQUANTIFIED_LABELS);
    }

    consensus.resolveUniqueIdConflicts();
    consensus.ensureUniqueId(); // the map's own id is reset when the fractions are combined
    addDataProcessing_(consensus, getProcessingInfo_(DataProcessing::QUANTITATION));

    //-------------------------------------------------------------
    // output
    //-------------------------------------------------------------
    if (!out_qpx.empty())
    {
      OPENMS_LOG_INFO << "Exporting QPX Parquet files to: " << out_qpx << endl;

      // every deterministic refusal of the three views, before the first file is written
      if (!QPXCollectionExport::requireExportable(consensus, design_))
      {
        return CANNOT_WRITE_OUTPUT_FILE; // already logged
      }
      // no partial collection is left behind without the commit below
      QPXCollectionExport::Transaction qpx_collection(out_qpx);

      QPXIdentity::FeatureLinks feature_links;
      if (!ConsensusMapArrowExport::exportToParquet(consensus, out_qpx + "/quantms.feature.parquet", ParquetWriteConfig{}, &feature_links))
      {
        OPENMS_LOG_ERROR << "Failed to write features Parquet file" << endl;
        return CANNOT_WRITE_OUTPUT_FILE;
      }

      PeptideIdentificationList all_pepids;
      for (const auto& feature : consensus)
      {
        for (const auto& pepid : feature.getPeptideIdentifications()) { all_pepids.push_back(pepid); }
      }
      for (const auto& pepid : consensus.getUnassignedPeptideIdentifications()) { all_pepids.push_back(pepid); }

      if (!QPXFile::exportToParquet(consensus.getProteinIdentifications(), all_pepids, out_qpx + "/quantms.psm.parquet",
                                    /*export_all_psms=*/false, ParquetWriteConfig{}, &feature_links))
      {
        OPENMS_LOG_ERROR << "Failed to write PSM Parquet file" << endl;
        return CANNOT_WRITE_OUTPUT_FILE;
      }

      if (!ProteinGroupArrowExport::exportToParquet(consensus, design_, out_qpx + "/quantms.pg.parquet"))
      {
        OPENMS_LOG_ERROR << "Failed to write protein groups Parquet file" << endl;
        return CANNOT_WRITE_OUTPUT_FILE;
      }

      qpx_collection.commit(); // all three views written
    }

    if (!out_cxml.empty())
    {
      FileHandler().storeConsensusFeatures(out_cxml, consensus, {FileTypes::CONSENSUSXML}, log_type_);
    }

    if (!out.empty())
    {
      const bool report_unidentified_features(false);
      const bool report_unmapped(true);
      const bool report_subfeatures(false);
      const bool report_unidentified_spectra(false);
      const bool report_not_only_best_psm_per_spectrum(false);
      MzTabFile().store(out, consensus, false, report_unidentified_features, report_unmapped, report_subfeatures,
                        report_unidentified_spectra, report_not_only_best_psm_per_spectrum);
    }

    return EXECUTION_OK;
  }

private:
  std::string labels_;                 ///< '-labels' as given
  int max_nr_labelled_aas_ = 0;        ///< '-max_nr_labelled_aas'
  Size multiplicity_ = 0;              ///< number of channels
  vector<std::string> channel_labels_; ///< column header label per channel, e.g. "no_label", "Lys8Arg10"
  vector<std::pair<std::string, std::string>> required_label_mods_; ///< (short, PSI-MS name) of every label the search must contain
  vector<set<std::string>> channel_label_sets_;  ///< per channel: its labels' short names ({"no_label"} for the light channel)
  vector<std::string> channel_descriptions_;     ///< per channel: PSI-MS names of its labels, for the column header
  std::map<std::string, vector<std::string>> label_long_to_shorts_; ///< PSI-MS name -> short names of '-labels' using it
  bool match_between_runs_ = false;    ///< keep unidentified multiplets so that linking can identify them from other runs
  bool knock_out_ = true;              ///< detect partial multiplets as well ('algorithm:knock_out', up to four channels)
  bool requantify_ = true;             ///< re-quantify the channels of partial multiplets from the MS1 peaks
  bool join_channels_ = true;          ///< join the separately detected channels of an identified peptide ('-join_channels')
  double requantify_min_correlation_ = 0.5; ///< minimum elution correlation of a re-quantified channel with its reference ('-requantify_min_correlation')
  double requantify_correlation_ratio_ = 0.3; ///< the correlation is checked below this fraction of the reference's intensity ('-requantify_correlation_ratio')
  bool targeted_quantification_ = true; ///< quantify identified peptides without a multiplet from their MS1 peaks ('-targeted_quantification')
  ExperimentalDesign design_;
  std::string picked_decoy_string_;
  bool picked_decoy_prefix_ = true;
  bool picked_decoy_seen_ = false;
};

int main(int argc, const char** argv)
{
  TOPPMS1LabeledWorkflow tool;
  return tool.main(argc, argv);
}

/// @endcond
