// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

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
#include <OpenMS/ANALYSIS/QUANTITATION/DDAWorkflowCommons.h>
#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
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
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/METADATA/MS1LabelState.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/SpectrumMetaDataLookup.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/SYSTEM/File.h>

#include "MS1LabeledRatioQuantifier.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <map>
#include <set>
#include <string>
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
  - Spectra in mzML format, one file per LC-MS run (@p in). Only MS1 spectra are used; profile and centroided
    data are both accepted (see @p algorithm:spectrum_type).
  - Identifications (@p ids), one file per spectra file in the same order, already filtered at PSM level
    (e.g. q-value < 0.01) and carrying Posterior Error Probability scores, e.g. produced with
    @ref TOPP_PercolatorAdapter (@p -score_type pep) or @ref TOPP_IDPosteriorErrorProbability, followed by
    @ref TOPP_IDFilter. Results from several search engines must be combined with @ref TOPP_ConsensusID
    rather than concatenated.
  - The label specification (@p labels), in the syntax of @ref TOPP_FeatureFinderMultiplex:
    <tt>[][Lys8,Arg10]</tt> for SILAC, <tt>[][Lys4,Arg6][Lys8,Arg10]</tt> for triple SILAC,
    <tt>[Dimethyl0][Dimethyl6]</tt> for Dimethyl. Every bracket is one channel; the channels are numbered
    from 1 in this order and are the @c Label column of the experimental design.
  - Optionally an experimental design (@p design) with the columns Fraction_Group, Fraction,
    Spectra_Filepath, Label and Sample (see @ref OpenMS::ExperimentalDesign "ExperimentalDesign"). One row per
    (file, channel). Without a design every file is an unfractionated fraction group and every
    (file, channel) is its own sample.
  - Optionally a protein database (@p fasta). The identifications are then re-indexed with
    @ref TOPP_PeptideIndexer, which annotates protein sequences (for coverage), the decoy status and the
    theoretical peptide uniqueness needed by <tt>-protein_quantification strictly_unique_peptides</tt>.

<b>Important:</b> the labels have to be part of the database search as (variable) modifications, e.g.
<tt>Label:13C(6)15N(2) (K)</tt> and <tt>Label:13C(6)15N(4) (R)</tt> for Lys8/Arg10. Otherwise the MS2 spectra of the
labeled channels stay unidentified and the observed mass shifts cannot be reconciled with the peptide
sequences, so most multiplets end up as conflicts. The tool checks the search parameters recorded in @p ids for
the modifications implied by @p labels and refuses to run if they are missing (use @p -force to proceed anyway).

<b>Workflow</b>
  -# Per run: detection of peptide multiplets in the MS1 data (@ref OpenMS::FeatureFinderMultiplexAlgorithm,
     parameters in @p algorithm and @p label_mass_shifts), mapping of the identifications onto the multiplets
     (@p id_mapping), reduction to one identification per multiplet, and consolidation of quantitative and
     sequence information (@ref OpenMS::MultiplexResolverAlgorithm, @p resolver): multiplets whose mass shifts
     contradict the labels found in the annotated sequence are removed from quantification (their
     identifications are kept for protein inference), incomplete multiplets are completed with dummy
     features (intensity 0 = absent, NaN = not quantifiable). Multiplets without identification are dropped,
     unless @p match_between_runs is set (see below). Once the resolver has used the labels, every identification
     is reduced to the <em>peptide identity</em>: the label modifications of @p labels are removed from the sequence,
     because the label belongs to the channel, not to the peptide, and the light and the heavy spectrum of one
     peptide have to name one peptide for linking, match between runs, inference and quantification (the
     convention MaxQuant uses as well). The label state stays documented on every identification as the meta
     values @c labeled_sequence (the peptidoform as searched), @c removed_labels (e.g. @c Lys8, or @c none) and
     @c label_channel (the 1-based channel the spectrum belongs to, i.e. the @c Label of the experimental
     design). The values also sit on every quantified consensus feature, for the identification it is quantified
     under. mzTab reports them as <tt>opt_global_*</tt> columns of the peptide and PSM sections, the QPX feature
     and psm views as @c cv_params. PSM-level output describes the spectrum match and therefore reports the
     peptidoform as searched (mzTab PSM section, QPX psm view), feature-level output the peptide identity
     (see @ref OpenMS::MS1LabelState). The column headers describe every channel's labels in @c channel_description.
     A spectrum match that was mapped onto several multiplets stays on the one whose matched channel is closest to
     the precursor; distinct spectra of one peptide on distinct multiplets are all kept, and their channel values
     add up per peptide and charge in the quantification.
  -# Per fraction: retention time alignment of the runs (@p alignment, identification-based, aligned to the
     run with most identifications) and linking of the multiplets across runs (@p linking); the channels of
     every run are kept as sub-features, so the linked map has one column per (run, channel). Fractions are
     linked separately and then combined column-wise, exactly like ProteomicsLFQ does; a fraction measured in a
     single run is passed through. With @p match_between_runs, unidentified multiplets take part in the linking
     and take over the identification of a multiplet at the same position in another run (the SILAC equivalent
     of ProteomicsLFQ's <tt>-targeted_only false</tt>); multiplets that stay unidentified are not quantified.
  -# Protein inference over all runs (@p protein_inference), protein (and optionally PSM/peptide) FDR
     filtering, and peptide and protein quantification (@p ProteinQuantification), where the fractions of a
     fraction group are aggregated according to the design.

<b>Ratios.</b> A labeled experiment measures its channels in one run, so its quantity is their <em>ratio</em>,
and the tool computes it the way MaxQuant does, as a median of ratios rather than as a ratio of aggregated
intensities (@p ratios, see @ref MS1LabeledRatioQuantifier):
  - per multiplet and run, the intensity of a channel over that of the reference channel (@p ratios:reference_channel,
    the light one by default). Only positive, finite channels take part: an absent (dummy, intensity 0) or
    not-quantifiable channel is no measurement of a ratio.
  - per peptide and fraction group, the median of its evidence ratios.
  - per protein group and fraction group, the median of its peptides' ratios, reported only from
    @p ratios:min_ratio_count peptides upwards (MaxQuant's "min. ratio count", 2 by default), with the number of
    contributing peptides next to it. Every ratio is also reported divided by the median peptide ratio of its
    (fraction group, channel), i.e. normalized on the assumption that most peptides do not change.

The reference channel is reported with the ratio 1.0 it has by construction, wherever another channel was
measured against it, so that every annotation covers the complete set of channels.

The ratios are annotated on the consensus features (@c evidence_ratio*, @c peptide_ratio*) and on the protein
groups, so they reach the consensusXML and the mzTab peptide section (as <tt>opt_global_*</tt> columns). In the
QPX pg view, whose rows are one per (protein group, fraction group, channel), they are written as that row's
@c additional_intensities, named @c ratio and @c ratio_normalized under the row's own channel label; the number
of contributing peptides sits in @c cv_params as @c ratio_count, being a count rather than an intensity. No
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
      "Input: spectra files (mzML), one per LC-MS run. Only MS1 spectra are used; profile and centroided data are accepted.");
    setValidFormats_("in", {"mzML"});

    registerInputFileList_("ids", "<file list>", StringList(),
      "Identifications filtered at PSM level (e.g., q-value < 0.01), one per spectra file in the same order.\n"
      "The identifications must carry Posterior Error Probability scores (e.g. PercolatorAdapter with -score_type pep,\n"
      "or IDPosteriorErrorProbability) and the labels must have been searched as (variable) modifications.\n"
      "Combine results from several search engines with ConsensusID rather than concatenating them.");
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
      "Cannot be combined with 'algorithm:knock_out': the channel order of an unidentified multiplet is only known from its detection pattern.", false);
    setValidStrings_("match_between_runs", {"true", "false"});

    // feature detection: the FeatureFinderMultiplex parameters, minus the two that are set at tool level
    Param ffm_defaults = FeatureFinderMultiplexAlgorithm().getDefaults();
    Param detection_defaults = ffm_defaults.copy("algorithm:", true);
    detection_defaults.remove("labels");
    detection_defaults.remove("max_nr_labelled_aas");
    Param mass_shift_defaults = ffm_defaults.copy("labels:", true);

    // resolver: tolerances only, labels and max_nr_labelled_aas are shared with the feature detection
    Param resolver_defaults = MultiplexResolverAlgorithm().getDefaults().copy("algorithm:", true);
    resolver_defaults.remove("labels");
    resolver_defaults.remove("max_nr_labelled_aas");

    Param idm_defaults = IDMapper().getDefaults();
    idm_defaults.addTag("mz_reference", "advanced");
    idm_defaults.addTag("ignore_charge", "advanced");

    Param ma_defaults = MapAlignmentAlgorithmIdentification().getDefaults();
    ma_defaults.setValue("max_rt_shift", 0.1);
    ma_defaults.setValue("use_unassigned_peptides", "false");
    ma_defaults.setValue("use_feature_rt", "true");
    for (const auto& s : {"score_type", "score_cutoff", "min_score", "use_unassigned_peptides", "use_feature_rt", "use_adducts"})
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
      OPENMS_LOG_WARN << "Warning: No mzML origin annotated in ID file. This can lead to errors or unexpected behaviour later: "
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
                         "measurement more than once. To control this upstream, combine search engine results with ConsensusID "
                         "(-algorithm best -keep_old_scores) rather than concatenating them." << endl;
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
  ExitCodes loadAndCleanupIDFile_(
    const std::string& id_file_abs_path,
    const std::string& mz_file,
    const std::string& in_db,
    const Size fraction_group,
    const Size fraction,
    vector<ProteinIdentification>& protein_ids,
    PeptideIdentificationList& peptide_ids,
    set<std::string>& fixed_modifications,    // adds to
    set<std::string>& variable_modifications) // adds to
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

    // make the run identifier unique across input files (users split mzML and id files before running the analysis)
    const std::string old_identifier = protein_ids[0].getIdentifier();
    const std::string new_identifier = old_identifier + "_" + StringUtils::toStr(fraction_group) + "F" + StringUtils::toStr(fraction);
    protein_ids[0].setIdentifier(new_identifier);
    for (PeptideIdentification& p : peptide_ids)
    {
      if (p.getIdentifier() == old_identifier) { p.setIdentifier(new_identifier); }
      else { OPENMS_LOG_WARN << "Peptide ID identifier found not present in the protein ID" << endl; }
    }

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
      SpectrumMetaDataLookup::addMissingSpectrumReferences(peptide_ids, mz_file_abs_path, true);
    }

    // deliberately last, i.e. after the spectrum-reference repair above
    reduceToOnePerSpectrum_(peptide_ids, id_file_abs_path);

    return EXECUTION_OK;
  }

  //-------------------------------------------------------------
  // per run: multiplet detection, ID mapping, resolving
  //-------------------------------------------------------------

  /// Detect the peptide multiplets of one run in its MS1 spectra (FeatureFinderMultiplex, incl. its FAIMS handling)
  ExitCodes detectMultiplets_(const std::string& mz_file, ConsensusMap& multiplets, MSExperiment& blacklist)
  {
    MSExperiment exp;
    FileHandler fh;
    fh.getOptions().setMSLevels({1}); // only MS1 spectra are needed
    fh.loadExperiment(mz_file, exp, {FileTypes::MZML}, log_type_);

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
      for (auto& spec : algorithm.getBlacklist())
      {
        if (has_faims) { spec.setMetaValue(Constants::UserParam::FAIMS_CV, group_cv); }
        blacklist.addSpectrum(std::move(spec));
      }
    }
    multiplets.ensureUniqueId();
    multiplets.sortByPosition();
    // the resolver looks peaks up by RT range, which needs the spectra of all FAIMS groups in one order
    blacklist.sortSpectra();
    blacklist.updateRanges();

    OPENMS_LOG_INFO << "Detected " << multiplets.size() << " peptide multiplet(s) in " << File::basename(mz_file) << "." << endl;
    return EXECUTION_OK;
  }

  /**
    @brief Attach every spectrum match to the closest channel, and to the closest multiplet only.

    IDMapper attaches an identification to every multiplet within its tolerances and records the
    first channel within tolerance as its 'map_index'. The resolver reads that index as the channel
    the spectrum belongs to, so it is set to the channel closest to the precursor m/z here (the same
    channel whenever the tolerance is smaller than the label shift). One spectrum is evidence for one
    multiplet, too (and a QPX PSM may reference one feature): it stays with the multiplet whose
    closest channel is closest to the precursor in m/z, then in RT. This runs before the best-scoring
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

        // the channel closest to the precursor is the one the spectrum belongs to
        double best_ppm = std::numeric_limits<double>::infinity(), best_rt = 0.0;
        Size best_map_index = 0;
        for (const FeatureHandle& handle : cf.getFeatures())
        {
          const double ppm = std::abs(handle.getMZ() - id.getMZ()) / handle.getMZ() * 1e6;
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
    ExitCodes e = loadAndCleanupIDFile_(id_file_abs_path, mz_file, in_db, fraction_group, fraction, protein_ids, peptide_ids,
                                        fixed_modifications, variable_modifications);
    if (e != EXECUTION_OK) { return e; }

    ConsensusMap multiplets;
    MSExperiment blacklist;
    e = detectMultiplets_(mz_file, multiplets, blacklist);
    if (e != EXECUTION_OK) { return e; }

    // map the identifications onto the multiplets: a multiplet matches if any of its channels does,
    // and the matched channel is recorded (needed by the resolver)
    {
      IDMapper mapper;
      Param idm_param = getParam_().copy("id_mapping:", true);
      writeDebug_("Parameters passed to IDMapper", idm_param, 3);
      mapper.setParameters(idm_param);
      mapper.annotate(multiplets, peptide_ids, protein_ids, /*measure_from_subelements=*/true, /*annotate_ids_with_subelements=*/true);
    }

    // a spectrum match that landed on several multiplets stays on the closest one only
    keepMultiAssignedOnClosest_(multiplets);

    // one identification (the best) per multiplet
    IDConflictResolverAlgorithm::resolve(multiplets);

    Size identified = 0;
    for (const auto& cf : multiplets) { identified += !cf.getPeptideIdentifications().empty(); }
    OPENMS_LOG_INFO << "Mapped identifications onto " << identified << " of " << multiplets.size() << " multiplet(s)." << endl;

    // remove quant/ID conflicts and complete the multiplets
    {
      MultiplexResolverAlgorithm resolver;
      Param res_param;
      res_param.insert("algorithm:", getParam_().copy("resolver:", true));
      res_param.setValue("algorithm:labels", labels_);
      res_param.setValue("algorithm:max_nr_labelled_aas", max_nr_labelled_aas_);
      res_param.insert("labels:", getParam_().copy("label_mass_shifts:", true));
      writeDebug_("Parameters passed to MultiplexResolverAlgorithm", res_param, 3);
      resolver.setParameters(res_param);

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
          // identification of a multiplet at the same position in another run. Its channels are in
          // pattern order, which is what the columns mean (guaranteed by 'knock_out' = false).
          if (match_between_runs_ && cf.size() == multiplicity_)
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
                        << "with the expected position) and are reported as missing." << endl;
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
      // the reference is determined from the data (the run with most identifications)
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
  void link_(vector<ConsensusMap>& maps, double max_alignment_diff, ConsensusMap& consensus_fraction)
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
    linker.group(maps, consensus_fraction);
    // the columns of the result are the channels of the input maps, not the input maps themselves
    linker.transferSubelements(maps, consensus_fraction);

    consensus_fraction.applyMemberFunction(&UniqueIdInterface::setUniqueId);
    consensus_fraction.sortPeptideIdentificationsByMapIndex();
    OPENMS_LOG_INFO << "Linked " << maps.size() << " run(s) into " << consensus_fraction.size() << " consensus multiplet(s)." << endl;
  }

  void alignAndLink_(vector<ConsensusMap>& maps, ConsensusMap& consensus_fraction)
  {
    if (maps.size() > 1)
    {
      vector<TransformationDescription> transformations;
      const double max_alignment_diff = align_(maps, transformations);
      for (Size i = 0; i < maps.size(); ++i)
      {
        try
        {
          MapAlignmentTransformer::transformRetentionTimes(maps[i], transformations[i]);
        }
        catch (Exception::IllegalArgument& e)
        {
          OPENMS_LOG_WARN << e.what() << endl;
        }
      }
      link_(maps, max_alignment_diff, consensus_fraction);
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

  /**
    @brief Reduce the hit's sequence to the peptide identity and record its label state.

    The label modifications of '-labels' are removed from the sequence (the channel, not the
    peptide, carries the label), and the hit is annotated with 'labeled_sequence' (as searched),
    'removed_labels' (short names, or 'none') and 'label_channel' (1-based; the 'Label' of the
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
    hit.setMetaValue(MS1LabelState::LABEL_CHANNEL, channelOfLabels_(removed));
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
    if (match_between_runs_ && getParam_().getValue("algorithm:knock_out").toString() == "true")
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "'-match_between_runs true' cannot be combined with '-algorithm:knock_out': with knock-outs, the channels of a multiplet "
        "are only ordered by its sequence, which an unidentified multiplet does not have.");
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

    // the design's channels must be the ones of '-labels'
    if (design_.getNumberOfLabels() != multiplicity_)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "The experimental design uses " + StringUtils::toStr(design_.getNumberOfLabels()) + " label(s), but '-labels " + labels_ + "' describes "
        + StringUtils::toStr(multiplicity_) + " channels. The 'Label' column has to number the channels of '-labels' from 1.");
    }

    // input files must be covered by the design, with every channel
    set<std::string> in_basenames;
    for (const auto& f : in) { in_basenames.insert(File::basename(f)); }
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
                      << " files from experimental design were not passed as mzMLs. Continuing with subset if the fractions still match." << endl;
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

    // mzML file -> id file (same order)
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
    quantifier.readQuantData(consensus, design_);
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
