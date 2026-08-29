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
#include <OpenMS/ANALYSIS/ID/IDBoostGraph.h>
#include <OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmTreeGuided.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PipEchoAlgorithm.h>
#include <OpenMS/ANALYSIS/QUANTITATION/DDAWorkflowCommons.h>
#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/CalibrationData.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>
#include <OpenMS/FEATUREFINDER/FeatureFindingMetabo.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/MSstatsFile.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ConversionHelper.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/SpectrumMetaDataLookup.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/FEATUREFINDER/Biosaur2Algorithm.h>
#ifdef WITH_OPENTIMS
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#endif
#include <OpenMS/ML/SVM/SimpleSVM.h>
#include <OpenMS/PROCESSING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/PROCESSING/CALIBRATION/MZTrafoModel.h>
#include <OpenMS/PROCESSING/CALIBRATION/PrecursorCorrection.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/FORMAT/ConsensusMapArrowExport.h>
#include <OpenMS/FORMAT/ProteinGroupArrowExport.h>
#include <OpenMS/FORMAT/QPXCollectionExport.h>
#include <OpenMS/FORMAT/QPXIdentity.h>
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/FORMAT/FeatureMapArrowIO.h>
#include <OpenMS/CONCEPT/VersionInfo.h>

#include <filesystem>
#include <iomanip>
#include <sstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_ProteomicsLFQ ProteomicsLFQ

ProteomicsLFQ performs label-free quantification of peptides and proteins. @n

Input: @n
  - Spectra in mzML format or Bruker .d directories (TimsTOF PASEF)
  - Identifications in idXML or mzIdentML format with posterior error probabilities
   as score type.
   To generate those we suggest to run:
    1. PeptideIndexer to annotate target and decoy information.
    2. PSMFeatureExtractor to annotate percolator features.
    3. PercolatorAdapter tool (score_type = 'q-value', -post-processing-tdc)
    4. IDFilter (pep:score = 0.01) to filter PSMs at 1% FDR

   Exactly one identification run per ID file is required, and merged ID runs are not supported.
   One identification per spectrum is expected as well: ProteomicsLFQ measures one value per
   (spectrum, peptidoform, charge), so where several identifications of one spectrum agree on
   all three, only the best-scoring one is kept and the reduction is reported. The others would
   otherwise count the same measurement more than once, in the PSM-level FDR and in every output.
   Results from several search engines must therefore be combined - with @ref TOPP_ConsensusID
   (@p -algorithm best @p -keep_old_scores, which preserves each engine's score) - rather than
   simply concatenated. Identifications of one spectrum that name *different* peptidoforms are
   left alone: a chimeric spectrum yields two distinct measurements.

  - An experimental design file: @n
   (see @ref OpenMS::ExperimentalDesign "ExperimentalDesign" for details) @n
  - A protein database in with appended decoy sequences in FASTA format @n
   (e.g., generated by the OpenMS DecoyDatabase tool) @n

Processing: @n
ProteomicsLFQ has different methods to extract features: ID-based (targeted only), or both ID-based and untargeted.
  1. The first method uses targeted feature dectection using RT and m/z information derived from identification data to extract features.
     Note: only identifications found in a particular MS run are used to extract features in the same run.
     No transfer of IDs (match between runs) is performed.
  2. The second method adds untargeted feature detection to obtain quantities from unidentified features.
     Transfer of Ids (match between runs) is performed by transfering feature identifications to coeluting, unidentified features with similar mass
and RT in other runs.

@b Resuming @b and @b distributing @b feature @b detection (@p -feat_dir): @n
Feature detection is the expensive part of the workflow and each MS run is detected independently
of every other; alignment, linking, inference and quantification need all runs at once. @p -feat_dir
names a directory of per-run feature checkpoints and applies one rule to every row of the
experimental design: reuse its checkpoint if a valid one is there, otherwise detect the run from
@p -in / @p -ids and write one, otherwise fail.

Everything follows from that rule:
@code
  // one machine, resumable -- interrupt it, run it again, it continues
  ProteomicsLFQ -design d.tsv -in *.mzML -ids *.idXML -fasta db.fasta -feat_dir ckpt/ -out r.mzTab

  // many machines: one detect-only invocation per run, in any order, no coordination
  ProteomicsLFQ -design d.tsv -in a.mzML -ids a.idXML -fasta db.fasta -feat_dir ckpt/ -detect_only

  // then combine, reading no mzML, no idXML and no run-level FASTA content at all
  ProteomicsLFQ -design d.tsv -fasta db.fasta -feat_dir ckpt/ -out r.mzTab
@endcode

A checkpoint records the parameters, OpenMS build, experimental-design row and input files it was
produced from, and is refused if any of those disagree with the run trying to use it - naming the
setting that differs. This is what makes reuse safe rather than merely convenient: nothing else
would stop half a study being detected with one setting and half with another. There is no way to
combine checkpoints that disagree: @p -force_recompute detects the affected runs again and rewrites
their checkpoints.

@p -feat_dir requires an explicit @p -design (a generated one would label every separately detected
run as the first) and @p -fasta (a checkpoint has to carry the peptide-indexing results, which a
combining run cannot reconstruct), and does not apply to @p spectral_counting.

Note on scale: the combining step holds every run's features of a fraction in memory at once, so its
ceiling is set by features per run rather than by run count. Measured on a 7-run dataset averaging
about 10,000 features per run, the marginal cost is roughly 2 kB per feature plus 10 MB per run, so
200 such runs in one fraction need about 9 GB; a deep-proteome experiment at ~60,000 features per run
would need several times that.

@b FAIMS (Field Asymmetric Ion Mobility Spectrometry): @n
FAIMS data is automatically detected based on compensation voltage (CV) annotations in the mzML file.
The data is split by CV and processed separately for each voltage group during feature detection.
Features representing the same analyte detected at different CV values are merged automatically.
The merged features are then aligned and linked across runs based on RT and m/z.
No special preparation of the input mzML file is required.

@b Bruker .d (TimsTOF PASEF): @n
Bruker .d directories containing DDA-PASEF data are supported directly.
When .d input is detected, the tool automatically:
  - Skips centroiding (PeakPickerHiRes) to preserve per-peak ion mobility data
  - Skips precursor mass correction (not IM-aware)
  - Forces Biosaur2Algorithm for seed generation (FeatureFinderMultiplex does not support IM_PEAK)
  - Estimates chromatographic FWHM from Biosaur2 feature extents
Identification files should be generated with SageAdapter, which annotates
ion mobility values in the idXML output. FeatureFinderIdentificationAlgorithm
uses these IM annotations for targeted 2D chromatogram extraction (m/z + IM windowing).
MS1 frames are IM-centroided during loading using BrukerTimsFile's built-in Sage algorithm,
collapsing ~245k raw peaks/frame into ~10k centroided peaks with summed intensity.
Biosaur2 defaults are tuned for ProteomicsLFQ (mini=500, minlh=3, pasefminlh=2).
On HeLa 50ng 5-min timsTOF gradient: 34k seeds, 80% model fit success, 2,809
peptides quantified (Spearman r=0.62 vs Sage LFQ), 75s runtime, 1.3 GB memory.

Normalization: @n
  - ProteomicsLFQ does NOT normalize. Every output - mzTab, consensusXML, QPX Parquet and MSstats -
    reports the same un-normalized abundances, whatever combination of output files is
    requested. Earlier versions applied median scaling to the consensus features unless
    -out_msstats was given, which made the meaning of an intensity depend on which
    other output file happened to be requested, and recorded nothing about the transform in any of
    them.
  - Choosing a normalization is a separate, explicit step. Three routes, in increasing distance
    from this tool: @n
    - @ref TOPP_ConsensusMapNormalizer on -out_cxml (median, quantile, robust regression or
      thresholded scaling); @n
    - the ProteinQuantification:consensus:normalize parameter, which scales peptide abundances so
      that the median of each (fraction group, label) assay matches the overall median - i.e. it
      normalizes at the assay level rather than per fraction; @n
    - the downstream tool itself. MSstats and comparable consumers normalize their input
      by design.

Output (at least one required; each output is optional individually):
  - mzTab file with analysis results (@p out)
  - MSstats file with analysis results for statistical downstream analysis in MSstats (@p out_msstats)
  - consensusXML file for visualization and further processing in OpenMS (@p out_cxml)
  - QPX Parquet collection (@p out_qpx)

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_ProteomicsLFQ.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_ProteomicsLFQ.html
 **/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class ProteomicsLFQ : public TOPPBase
{
public:
  ProteomicsLFQ(): TOPPBase("ProteomicsLFQ", "A standard proteomics LFQ pipeline.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<file list>", StringList(),
      "Input files. Optional only when '-feat_dir' supplies a checkpoint for every run of the design.",
      false);
    setValidFormats_("in", ListUtils::create<std::string>("mzML"
#ifdef WITH_OPENTIMS
      ",d"
#endif
#ifdef WITH_THERMO_RAW
      ",raw"
#endif
    ));
    registerInputFileList_("ids", "<file list>", StringList(),
      "Identifications filtered at PSM level (e.g., q-value < 0.01)."
      "And annotated with PEP as main score.\n"
      "We suggest using:\n"
      "1. PSMFeatureExtractor to annotate percolator features.\n"
      "2. PercolatorAdapter tool (score_type = 'q-value', -post-processing-tdc)\n"
      "3. IDFilter (pep:score = 0.05)\n"
      "To obtain well calibrated PEPs and an initial reduction of PSMs\n"
      "ID files must be provided in same order as spectra files.\n"
      "One identification per spectrum is expected: where several identifications of one\n"
      "spectrum agree on peptidoform and charge, only the best-scoring one is kept, since\n"
      "the others would count the same measurement more than once. Combine results from\n"
      "several search engines with ConsensusID rather than concatenating them.", false);
    setValidFormats_("ids", ListUtils::create<std::string>("idXML,mzId,idparquet"));

    registerInputFile_("design", "<file>", "", "design file", false);
    setValidFormats_("design", ListUtils::create<std::string>("tsv"));

    registerInputFile_("fasta", "<file>", "", "fasta file", false);
    setValidFormats_("fasta", ListUtils::create<std::string>("fasta"));

    registerOutputFile_("out", "<file>", "", "Optional output mzTab file. At least one output must be specified.", false, false);
    setValidFormats_("out", ListUtils::create<std::string>("mzTab"));

    registerOutputFile_("out_msstats", "<file>", "", "Optional output MSstats input file. At least one output must be specified.", false, false);
    setValidFormats_("out_msstats", ListUtils::create<std::string>("csv"));


    registerOutputFile_("out_cxml", "<file>", "", "Optional output consensusXML file. At least one output must be specified.", false, false);
    setValidFormats_("out_cxml", ListUtils::create<std::string>("consensusXML"));

    registerOutputDir_("out_qpx", "<directory>", "", "Optional output directory for QPX Parquet files (quantms.feature.parquet, quantms.psm.parquet, quantms.pg.parquet). At least one output must be specified.", false, false);

    registerOutputDir_("feat_dir", "<directory>", "",
      "Directory of per-run feature checkpoints. For every run of the experimental design, a valid "
      "checkpoint here is reused instead of detecting features again; a run without one is detected "
      "from '-in'/'-ids' and its checkpoint written. This makes a run resumable, and lets the per-run "
      "work be distributed: run with '-detect_only' on each machine, then once over the design with "
      "neither '-in' nor '-ids'. Requires '-design' and '-fasta'.", false, false);
    registerFlag_("detect_only",
      "Stop after the per-run feature checkpoints have been written. No alignment, linking, inference "
      "or quantification is performed, and no result file is required. Requires '-feat_dir'.", false);
    registerFlag_("force_recompute",
      "Ignore existing checkpoints in '-feat_dir' and detect every run again, overwriting them.", true);

    registerDoubleOption_("proteinFDR", "<threshold>", 0.05, "Protein FDR threshold (0.05=5%).", false);
    setMinFloat_("proteinFDR", 0.0);
    setMaxFloat_("proteinFDR", 1.0);

    // TODO test rigorously
    registerStringOption_("picked_proteinFDR", "<choice>", "false", "Use a picked protein FDR?", false);
    setValidStrings_("picked_proteinFDR", {"true", "false"});

    registerDoubleOption_(
      "psmFDR", "<threshold>", 1.0,
      "FDR threshold for sub-protein level (e.g. 0.05=5%). Use -FDR_type to choose the level. Cutoff is applied at the highest level."
      " If Bayesian inference was chosen, it is equivalent with a peptide FDR",
      false);
    setMinFloat_("psmFDR", 0.0);
    setMaxFloat_("psmFDR", 1.0);

    registerStringOption_("FDR_type", "<threshold>", "PSM", "Sub-protein FDR level. PSM, PSM+peptide (best PSM q-value).", false);
    setValidStrings_("FDR_type", {"PSM", "PSM+peptide"});

    // TODO expose all parameters of the inference algorithms (e.g. aggregation methods etc.)?
    registerStringOption_("protein_inference", "<option>", "aggregation",
      "Infer proteins:\n"
      "aggregation  = aggregates all peptide scores across a protein (using the best score) \n"
      "bayesian     = computes a posterior probability for every protein based on a Bayesian network.\n"
      "               Note: 'bayesian' only uses and reports the best PSM per peptide.",
      false, true);
    setValidStrings_("protein_inference", ListUtils::create<std::string>("aggregation,bayesian"));

    registerStringOption_("protein_quantification", "<option>", "unique_peptides",
      "Quantify proteins based on:\n"
      "unique_peptides = use peptides mapping to single proteins or a group of indistinguishable proteins"
      "(according to the set of experimentally identified peptides).\n"
      "strictly_unique_peptides = use peptides mapping to a unique single protein only.\n"
      "shared_peptides = use shared peptides only for its best group (by inference score)", false, true);
    setValidStrings_("protein_quantification", ListUtils::create<std::string>("unique_peptides,strictly_unique_peptides,shared_peptides"));
    registerStringOption_("quantification_method", "<option>",
      "feature_intensity",
      "feature_intensity: MS1 signal.\n"
      "spectral_counting: PSM counts.", false, false);
    setValidStrings_("quantification_method", ListUtils::create<std::string>("feature_intensity,spectral_counting"));

    registerStringOption_("targeted_only", "<option>", "false",
      "true: Only ID based quantification.\n"
      "false: include unidentified features so they can be linked to identified ones (=match between runs).", false, false);
    setValidStrings_("targeted_only", ListUtils::create<std::string>("true,false"));

    registerStringOption_("pip_echo", "<option>", "false", "Perform match between runs (MBR) via PIP-ECHO", false, false);
    setValidStrings_("pip_echo", ListUtils::create<std::string>("true,false"));

    registerDoubleOption_(
      "feature_with_id_min_score", "<p-value>", 0.0,
      "The minimum probability (e.g.: 0.25) an identified (=id targeted) feature must have to be kept for alignment and linking (0=no filter).",
      false, true);
    setMinFloat_("feature_with_id_min_score", 0.0);
    setMaxFloat_("feature_with_id_min_score", 1.0);

    registerDoubleOption_(
      "feature_without_id_min_score", "<p-value>", 0.0,
      "The minimum probability (e.g.: 0.75) an unidentified feature must have to be kept for alignment and linking (0=no filter).", false, true);
    setMinFloat_("feature_without_id_min_score", 0.0);
    setMaxFloat_("feature_without_id_min_score", 1.0);

    registerStringOption_("mass_recalibration", "<option>", "false", "Mass recalibration.", false, true);
    setValidStrings_("mass_recalibration", ListUtils::create<std::string>("true,false"));

    registerStringOption_("alignment_order", "<option>", "star", "If star, aligns all maps to the reference with most IDs. If tree_guided, aligns maps in tree order (most similar pairs first).", false, true);
    setValidStrings_("alignment_order", ListUtils::create<std::string>("star,tree_guided"));

    registerStringOption_("keep_feature_top_psm_only", "<option>", "true", "If false, also keeps lower ranked PSMs that have the top-scoring"
                                                                     " sequence as a candidate per feature in the same file.", false, true);
    setValidStrings_("keep_feature_top_psm_only", ListUtils::create<std::string>("true,false"));

    registerTOPPSubsection_("Seeding", "Parameters for seeding of untargeted features");
    registerDoubleOption_("Seeding:intThreshold", "<threshold>", 1e4, "Peak intensity threshold applied in seed detection.", false, true);
    registerStringOption_("Seeding:charge", "<minChg:maxChg>", "2:5", "Charge range considered for untargeted feature seeds.", false, true); //TODO infer from IDs?
    registerDoubleOption_("Seeding:traceRTTolerance", "<tolerance(sec)>", 3.0, "Combines all spectra in the tolerance window to stabilize identification of isotope patterns. Controls sensitivity (low value) vs. specificity (high value) of feature seeds.", false, true); //TODO infer from average MS1 cycle time?
    registerStringOption_("Seeding:algorithm", "<choice>", "multiplex",
      "Algorithm for untargeted seed feature detection.\n"
      "multiplex: FeatureFinderMultiplexAlgorithm (default, current behavior).\n"
      "biosaur2: Biosaur2Algorithm (handles IM_PEAK/PASEF data natively).",
      false, false);
    setValidStrings_("Seeding:algorithm", {"multiplex", "biosaur2"});

    /// TODO: think about export of quality control files (qcML?)

    Param pp_defaults = PeakPickerHiRes().getDefaults();
    for (const auto& s : {"report_FWHM", "report_FWHM_unit", "SignalToNoise:win_len", "SignalToNoise:bin_count",
                          "SignalToNoise:min_required_elements", "SignalToNoise:write_log_messages"})
    {
      pp_defaults.addTag(s, "advanced");
    }

    Param ffi_defaults = FeatureFinderIdentificationAlgorithm().getDefaults();
    ffi_defaults.setValue("svm:samples", 10000); // restrict number of samples for training
    ffi_defaults.setValue("svm:log2_C", DoubleList({-2.0, 5.0, 15.0}));
    ffi_defaults.setValue("svm:log2_gamma", DoubleList({-3.0, -1.0, 2.0}));
    ffi_defaults.setValue("svm:min_prob", 0.9); // keep only feature candidates with > 0.9 probability of correctness

    // hide entries
    for (const auto& s : {"svm:samples", "svm:log2_C", "svm:log2_gamma", "svm:min_prob", "svm:no_selection", "svm:xval_out", "svm:kernel", "svm:xval",
                          "candidates_out", "extract:n_isotopes", "model:type"})
    {
      ffi_defaults.addTag(s, "advanced");
    }
    ffi_defaults.remove("detect:peak_width"); // set from data

    Param ma_defaults = MapAlignmentAlgorithmTreeGuided().getDefaults();

    ma_defaults.setValue("align_algorithm:max_rt_shift", 0.1);
    ma_defaults.setValue("align_algorithm:use_unassigned_peptides", "false");
    ma_defaults.setValue("align_algorithm:use_feature_rt", "true");

    // hide entries
    for (const auto& s :
         {"align_algorithm:use_unassigned_peptides", "align_algorithm:use_feature_rt", "align_algorithm:score_cutoff", "align_algorithm:min_score"})
    {
      ma_defaults.addTag(s, "advanced");
    }

    // Param fl_defaults = FeatureGroupingAlgorithmKD().getDefaults();
    Param fl_defaults = FeatureGroupingAlgorithmQT().getDefaults();
    fl_defaults.setValue("distance_MZ:max_difference", 10.0);
    fl_defaults.setValue("distance_MZ:unit", "ppm");
    fl_defaults.setValue("distance_MZ:weight", 5.0);
    fl_defaults.setValue("distance_intensity:weight", 0.1);
    fl_defaults.setValue("use_identifications", "true");
    fl_defaults.remove("distance_RT:max_difference"); // estimated from data
    for (const auto& s : {"distance_MZ:weight", "distance_intensity:weight", "use_identifications", "ignore_charge", "ignore_adduct"})
    {
      fl_defaults.addTag(s, "advanced");
    }

    // For PIP-ECHO:
    Param pip_echo_defaults = PipEchoAlgorithm().getDefaults();
    pip_echo_defaults.remove("distance_RT:max_difference"); // estimated from data

    Param pq_defaults = PeptideAndProteinQuant().getDefaults();
    // overwrite algorithm default, so we export everything (important for copying back MSstats results)
    pq_defaults.setValue("top:include_all", "true");
    pq_defaults.addTag("top:include_all", "advanced");

    Param bio_defaults = Biosaur2Algorithm().getDefaults();
    bio_defaults.setValue("mini", 500.0);   // filter low-intensity noise peaks (default 1.0 too permissive)
    bio_defaults.setValue("minlh", 3);      // require hills spanning >= 3 scans (default 2 keeps transient noise)
    bio_defaults.setValue("pasefminlh", 2); // require >= 2 raw points per PASEF cluster (default 1)
    for (auto it = bio_defaults.begin(); it != bio_defaults.end(); ++it)
    {
      bio_defaults.addTag(it.getName(), "advanced");
    }

    // combine parameters of the individual algorithms
    Param combined;
    combined.insert("Centroiding:", pp_defaults);
    combined.insert("PeptideQuantification:", ffi_defaults);
    combined.insert("Alignment:", ma_defaults);
    combined.insert("Linking:", fl_defaults);
    combined.insert("ProteinQuantification:", pq_defaults);
    combined.insert("Seeding:Biosaur2:", bio_defaults);
    combined.insert("PipEcho:", pip_echo_defaults);

    registerFullParam_(combined);
  }

  ExitCodes loadAndPreprocess_(const std::string& mz_file, MSExperiment& ms_out, bool& is_im_peak_data)
  {
#ifdef WITH_OPENTIMS
    if (FileHandler::getType(mz_file) == FileTypes::BRUKER_TDF)
    {
      // .d path: load with built-in IM centroiding (Sage algorithm, Lazear 2023)
      // This collapses raw IM profiles into single centroided peaks with summed intensity,
      // matching what Sage v0.15 does internally before LFQ feature tracing.
      BrukerTimsFile tdf;
      tdf.setLogType(log_type_);
      BrukerTimsFile::Config config;
      config.ms1_centroid_mz_ppm = 5.0f;  // same as Sage default
      config.ms1_centroid_im_pct = 3.0f;  // same as Sage default
      tdf.load(mz_file, ms_out, config);
      ms_out.updateRanges();

      if (ms_out.empty())
      {
        OPENMS_LOG_WARN << "The given file does not contain any spectra.";
        return INCOMPATIBLE_INPUT_DATA;
      }

      // Remove MS2 spectra entirely and filter MS1 spectra without IM data.
      // MS2 peak data is not needed (IDs come from external search engine).
      // MS1 spectra without IM arrays cause ChromatogramExtractorAlgorithm to throw
      // during IM-windowed extraction.
      auto& spectra = ms_out.getSpectra();
      for (auto& spec : spectra)
      {
        if (!spec.isSorted())
        {
          spec.sortByPosition();
        }
      }
      spectra.erase(
        std::remove_if(spectra.begin(), spectra.end(),
          [](const MSSpectrum& spec)
          {
            return spec.getMSLevel() != 1 || !spec.containsIMData();
          }),
        spectra.end());

      // No PeakPickerHiRes — Bruker TOF data is centroid-like; IM arrays must be preserved
      // No PrecursorCorrection — findHighestInWindow() is IM-unaware
      // No clearMetaDataArrays — IM per-peak float arrays must survive
      is_im_peak_data = true;
      OPENMS_LOG_INFO << "Loaded Bruker .d file with IM_PEAK data. Skipping centroiding and precursor correction.\n";
      return EXECUTION_OK;
    }
    else
#endif
    {
      // mzML path: existing centroid + precursor correction
      is_im_peak_data = false;
      return centroidAndCorrectPrecursors_(mz_file, ms_out);
    }
  }

  ExitCodes centroidAndCorrectPrecursors_(const std::string & mz_file, MSExperiment & ms_centroided)
  {
    Param pp_param = getParam_().copy("Centroiding:", true);
    writeDebug_("Parameters passed to PeakPickerHiRes algorithm", pp_param, 3);

    // create scope for raw data, so it is properly freed (Note: clear() is not sufficient)
    // load raw file

    PeakMap ms_raw;
    FileHandler().loadExperiment(mz_file, ms_raw, {FileTypes::MZML, FileTypes::BRUKER_TDF, FileTypes::RAW}, log_type_);
    ms_raw.clearMetaDataArrays();
    ms_raw.updateRanges();

    if (ms_raw.empty())
    {
      OPENMS_LOG_WARN << "The given file does not contain any spectra.";
      return INCOMPATIBLE_INPUT_DATA;
    }

    // remove MS2 peak data and check if spectra are sorted
    // TODO can we load just MS1 or do we need precursor information?
    for (auto& spec : ms_raw)
    {
      if (spec.getMSLevel() == 2)
      {
        spec.clear(false); // delete MS2 peaks
      }
      if (! spec.isSorted())
      {
        spec.sortByPosition();
        writeLogInfo_("Info: Sorted peaks by m/z.");
      }
    }

    //-------------------------------------------------------------
    // Centroiding of MS1
    //-------------------------------------------------------------
    PeakPickerHiRes pp;
    pp.setLogType(log_type_);
    pp.setParameters(pp_param);
    pp.pickExperiment(ms_raw, ms_centroided, true);

    //-------------------------------------------------------------
    // HighRes Precursor Mass Correction
    //-------------------------------------------------------------
    std::vector<double> deltaMZs, mzs, rts;
    std::set<Size> corrected_to_highest_intensity_peak = PrecursorCorrection::correctToHighestIntensityMS1Peak(
      ms_centroided,
      0.01, // check if we can estimate this from data (here it is given in m/z not ppm)
      false, // is ppm = false
      deltaMZs,
      mzs,
      rts
      );
    writeLogInfo_("Info: Corrected " + StringUtils::toStr(corrected_to_highest_intensity_peak.size()) + " precursors.");
    if (!deltaMZs.empty())
    {
      vector<double> deltaMZs_ppm, deltaMZs_ppmabs;
      for (Size i = 0; i != deltaMZs.size(); ++i)
      {
        deltaMZs_ppm.push_back(Math::getPPM(mzs[i], mzs[i] + deltaMZs[i]));
        deltaMZs_ppmabs.push_back(Math::getPPMAbs(mzs[i], mzs[i] + deltaMZs[i]));
      }

      double median = Math::median(deltaMZs_ppm.begin(), deltaMZs_ppm.end());
      double MAD = Math::MAD(deltaMZs_ppm.begin(), deltaMZs_ppm.end(), median);
      double median_abs = Math::median(deltaMZs_ppmabs.begin(), deltaMZs_ppmabs.end());
      double MAD_abs = Math::MAD(deltaMZs_ppmabs.begin(), deltaMZs_ppmabs.end(), median_abs);
      writeLogInfo_("Precursor correction:\n  median        = "
        + StringUtils::toStr(median) + " ppm  MAD = " + StringUtils::toStr(MAD)
        + "\n  median (abs.) = " + StringUtils::toStr(median_abs)
        + " ppm  MAD = " + StringUtils::toStr(MAD_abs));
    }
    return EXECUTION_OK;
  }

  // aligns the feature maps
  double align_(vector<FeatureMap>& feature_maps, vector<TransformationDescription>& transformations)
  {
    if (feature_maps.size() > 1) // do we have several maps to align / link?
    {
      Param mat_param = getParam_().copy("Alignment:", true);
      writeDebug_("Parameters passed to MapAlignmentAlgorithms", mat_param, 3);

      Param model_params = MapAlignerBase::getModelDefaults("b_spline");
      std::string model_type = model_params.getValue("type").toString();
      model_params = model_params.copy(model_type + ":", true);

      try
      {
        if (getStringOption_("alignment_order") == "star")
        {
          // Determine reference from data, otherwise a change in order of input files
          // leads to slightly different results
          const int reference_index(-1); // set no reference (determine from data)
          Param ma_param = mat_param.copy("align_algorithm:", true);
          writeDebug_("Parameters passed to MapAlignerIdentification", ma_param, 3);
          MapAlignmentAlgorithmIdentification aligner;
          aligner.setLogType(log_type_);
          aligner.setParameters(ma_param);
          aligner.align(feature_maps, transformations, reference_index);
        }
        else // tree-guided
        {
          MapAlignmentAlgorithmTreeGuided aligner;
          aligner.setLogType(log_type_);
          aligner.setParameters(mat_param);
          aligner.align(feature_maps, transformations);
        }
      }
      catch (Exception::MissingInformation& err)
      {
        if (getFlag_("force"))
        {
          OPENMS_LOG_ERROR << "Error: alignment failed. Details:\n"
                           << err.what() << "\nProcessing will continue using 'identity' transformations." << endl;
          model_type = "identity";
          transformations.resize(feature_maps.size());
        }
        else
          throw;
      }

      // find model parameters (if model_type == "identity" the fit is a NOP):
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

      // determine maximum RT shift after transformation that includes all high confidence IDs
      using TrafoStat = TransformationDescription::TransformationStatistics;
      for (auto& s : alignment_stats)
      {
        OPENMS_LOG_INFO << "Alignment differences (second) for percentiles (before & after): " << endl;
        OPENMS_LOG_INFO << ListUtils::concatenate(s.percents, "%\t") << "%" << endl;
        OPENMS_LOG_INFO << "before alignment:" << endl;
        for (const auto& p : s.percents)
        {
          OPENMS_LOG_INFO << (int)s.percentiles_before[p] << "\t";
        }
        OPENMS_LOG_INFO << endl;

        OPENMS_LOG_INFO << "after alignment:" << endl;
        for (const auto& p : s.percents)
        {
          OPENMS_LOG_INFO << (int)s.percentiles_after[p] << "\t";
        }
        OPENMS_LOG_INFO << endl;
      }

      double max_alignment_diff = std::max_element(alignment_stats.begin(), alignment_stats.end(), [](TrafoStat a, TrafoStat b) {
                                    return a.percentiles_after[100] < b.percentiles_after[100];
                                  })->percentiles_after[100];
      // sometimes, very good alignments might lead to bad overall performance. Choose 2 minutes as minimum.
      OPENMS_LOG_INFO << "Max alignment difference (seconds): " << max_alignment_diff << endl;
      max_alignment_diff = std::max(max_alignment_diff, 120.0); // minimum 2 minutes
      max_alignment_diff = std::min(max_alignment_diff, 600.0); // maximum 10 minutes
      return max_alignment_diff;
    }
    return 0;
  }

  void transform_(vector<FeatureMap>& feature_maps, vector<TransformationDescription>& transformations)
  {
    if (feature_maps.size() > 1 && ! transformations.empty())
    {
      // Apply transformations
      for (Size i = 0; i < feature_maps.size(); ++i)
      {
        try
        {
          MapAlignmentTransformer::transformRetentionTimes(feature_maps[i], transformations[i]);
        }
        catch (Exception::IllegalArgument& e)
        {
          OPENMS_LOG_WARN << e.what() << endl;
        }

        if (debug_level_ > 666)
        {
          // plot with e.g.:
          // Rscript ../share/OpenMS/SCRIPTS/plot_trafo.R debug_trafo_1.trafoXML debug_trafo_1.pdf
          FileHandler().storeTransformations("debug_trafo_" + StringUtils::toStr(i) + ".trafoXML", transformations[i], {FileTypes::TRANSFORMATIONXML});
        }
      }
    }
  }

  /**
   * Link all features of the given fraction.
   *
   * In other words, link multiple runs together into a ConsensusMap
   * which, if seeding was performed earlier, will lead to ID transfer
   * to those runs that are missing MS2 peaks.
   */
  void link_(vector<FeatureMap>& feature_maps, double median_fwhm, double max_alignment_diff, ConsensusMap& consensus_fraction)
  {
    //since requantification only happens with 2+ maps, we do not need to check/skip,
    //in case of a singleton fraction. Would throw an exception in linker.group

    Param fl_param = getParam_().copy("Linking:", true);
    writeDebug_("Parameters passed to feature grouping algorithm", fl_param, 3);

    writeDebug_("Linking: " + StringUtils::toStr(feature_maps.size()) + " features.", 1);

    if (getStringOption_("pip_echo") != "false")
    {
      PipEchoAlgorithm linker;

      Param pe_param = getParam_().copy("PipEcho:", true);
      pe_param.setValue("distance_RT:max_difference", 2.0 * max_alignment_diff + 2.0 * median_fwhm);
      // Pass the host's raw-spectra chromatographic FWHM down for PIP-ECHO's
      // local_rt:auto window sizing (used only when local_rt:enabled+auto are set).
      pe_param.setValue("local_rt:median_fwhm", median_fwhm);
      writeDebug_("Parameters passed to the PIP-ECHO algorithm", pe_param, 3);

      linker.setParameters(pe_param);
      linker.group(feature_maps, consensus_fraction);
    }
    else
    {
      FeatureGroupingAlgorithmQT linker;

      Param fl_param = getParam_().copy("Linking:", true);
      fl_param.setValue("distance_RT:max_difference", 2.0 * max_alignment_diff + 2.0 * median_fwhm);
      writeDebug_("Parameters passed to feature grouping algorithm", fl_param, 3);

      linker.setParameters(fl_param);
      linker.group(feature_maps, consensus_fraction);
    }

    OPENMS_LOG_INFO << "Size of consensus fraction: " << consensus_fraction.size() << endl;
    assert(! consensus_fraction.empty());
  }

  /// Align and link.
  void alignAndLink_(vector<FeatureMap>& feature_maps,
                     ConsensusMap& consensus_fraction,
                     vector<TransformationDescription>& transformations,
                     const double median_fwhm)
  {
    double max_alignment_diff(0.0);

    if (feature_maps.size() > 1)
    {
      max_alignment_diff = align_(feature_maps, transformations);

      transform_(feature_maps, transformations);

      link_(feature_maps, median_fwhm, max_alignment_diff, consensus_fraction);
    }
    else // only one feature map
    {
      MapConversion::convert(0, feature_maps.back(), consensus_fraction);
    }
  }

  /// determine in which runs of the current fraction a peptide was quantified
  /// returns map sequence+charge -> map index in consensus map that have non-zero quant values
  map<pair<std::string, UInt>, vector<int> > getPeptideOccurrence_(const ConsensusMap &cons)
  {
    map<Size, UInt> num_consfeat_of_size;
    map<Size, UInt> num_consfeat_of_size_with_id;

    map<pair<std::string, UInt>, vector<int> > seq_charge2map_occurence;
    for (const ConsensusFeature& cfeature : cons)
    {
      ++num_consfeat_of_size[cfeature.size()];
      const auto& pids = cfeature.getPeptideIdentifications();
      if (! pids.empty())
      {
        ++num_consfeat_of_size_with_id[cfeature.size()];

        // count how often a peptide/charge pair has been observed in the different maps
        const vector<PeptideHit>& phits = pids[0].getHits();
        if (! phits.empty())
        {
          const std::string s = phits[0].getSequence().toString();
          const int z = phits[0].getCharge();

          if (seq_charge2map_occurence[make_pair(s, z)].empty())
          {
            seq_charge2map_occurence[make_pair(s, z)] = vector<int>(cons.getColumnHeaders().size(), 0);
          }

          // assign id to all dimensions in the consensus feature
          for (auto const& f : cfeature.getFeatures())
          {
            Size map_index = f.getMapIndex();
            seq_charge2map_occurence[make_pair(s, z)][map_index] += 1;
          }
        }
      }
    }
    return seq_charge2map_occurence;
  }

  ExitCodes checkSingleRunPerID_(const vector<ProteinIdentification>& protein_ids, const std::string& id_file_abs_path)
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
      OPENMS_LOG_FATAL_ERROR << "ProteomicsLFQ does not support merged ID runs. ID file: " << id_file_abs_path << endl;
      return ExitCodes::INCOMPATIBLE_INPUT_DATA;
    }
    if (run_paths.empty())
    {
      OPENMS_LOG_WARN << "Warning: No mzML origin annotated in ID file. This can lead to errors or unexpected behaviour later: " << id_file_abs_path
                      << endl;
    }

    return EXECUTION_OK;
  }

  ExitCodes switchScoreType_(PeptideIdentificationList& peptide_ids, const std::string& id_file_abs_path)
  {
    // Check if score types are valid. TODO
    try
    {
      IDScoreSwitcherAlgorithm switcher;
      Size c = 0;
      switcher.switchToGeneralScoreType(peptide_ids, IDScoreSwitcherAlgorithm::ScoreType::PEP, c);
    }
    catch (Exception::MissingInformation&)
    {
      OPENMS_LOG_FATAL_ERROR << "ProteomicsLFQ expects a Posterior Error Probability score in all Peptide IDs. ID file: " << id_file_abs_path << endl;
      return ExitCodes::INCOMPATIBLE_INPUT_DATA;
    }
    return EXECUTION_OK;
  }

  /**
    @brief Reconcile one input file's decoy affix into the single study-wide one.

    Picked protein FDR takes one affix for the whole study, but it is inferred once per input
    file. This used to be a plain overwrite, so the last file silently decided it. The affix is
    a property of the FASTA and so is identical for every file of a well-formed run; two files
    disagreeing means one was searched against a different database, and the picked FDR would
    then be computed against an affix that does not match half the accessions. Keep the first
    and say so.
  */
  void recordDecoyAffix_(const std::string& decoy_string, bool decoy_prefix, bool decoy_inferred,
                         const std::string& source_file)
  {
    if (! decoy_inferred) return; // no FASTA given, nothing was inferred from this file

    if (! picked_decoy_seen_)
    {
      picked_decoy_string_ = decoy_string;
      picked_decoy_prefix_ = decoy_prefix;
      picked_decoy_seen_ = true;
    }
    else if (decoy_string != picked_decoy_string_ || decoy_prefix != picked_decoy_prefix_)
    {
      OPENMS_LOG_WARN << "Warning: decoy affix inferred from " << source_file << " ('"
                      << decoy_string << "', " << (decoy_prefix ? "prefix" : "suffix")
                      << ") differs from the one inferred from earlier input ('" << picked_decoy_string_
                      << "', " << (picked_decoy_prefix_ ? "prefix" : "suffix")
                      << "). Keeping the first. Picked protein FDR assumes one affix for the whole "
                         "study - check that all inputs were searched against the same database.\n";
    }
  }

  ExitCodes loadAndCleanupIDFile_(
    const std::string& id_file_abs_path,
    const std::string& mz_file,
    const std::string& in_db,
    const Size& fraction_group,
    const Size& fraction,
    vector<ProteinIdentification>& protein_ids,
    PeptideIdentificationList& peptide_ids,
    set<std::string>& fixed_modifications,  // adds to
    set<std::string>& variable_modifications, // adds to
    std::string& out_decoy_string,   // decoy affix this file's IDs imply
    bool& out_decoy_prefix,
    bool& out_decoy_inferred)        // false when no FASTA was given, so nothing was inferred
  {

    const std::string& mz_file_abs_path = File::absolutePath(mz_file);
    FileHandler().loadIdentifications(id_file_abs_path, protein_ids, peptide_ids,
        {FileTypes::IDXML, FileTypes::MZIDENTML, FileTypes::IDPARQUET}, log_type_);

    ExitCodes e = checkSingleRunPerID_(protein_ids, id_file_abs_path);
    if (e != EXECUTION_OK) return e;

    // Re-index
    if (! in_db.empty())
    {
      PeptideIndexing indexer;
      Param param_pi = indexer.getParameters();
      param_pi.setValue("missing_decoy_action", "silent");
      param_pi.setValue("write_protein_sequence", "true");
      param_pi.setValue("write_protein_description", "true");
      indexer.setParameters(param_pi);

      // stream data in fasta file
      FASTAContainer<TFI_File> fasta_db(in_db);
      PeptideIndexing::ExitCodes indexer_exit = indexer.run(fasta_db, protein_ids, peptide_ids);

      // The affix this file's identifications imply. It is reported back rather than recorded
      // here: picked protein FDR needs one affix for the whole study, and reconciling per-file
      // values into that one is the caller's business (see recordDecoyAffix_).
      out_decoy_string = indexer.getDecoyString();
      out_decoy_prefix = indexer.isPrefix();
      out_decoy_inferred = true;
      if ((indexer_exit != PeptideIndexing::ExitCodes::EXECUTION_OK) && (indexer_exit != PeptideIndexing::ExitCodes::PEPTIDE_IDS_EMPTY))
      {
        if (indexer_exit == PeptideIndexing::ExitCodes::DATABASE_EMPTY) { return INPUT_FILE_EMPTY; }
        else if (indexer_exit == PeptideIndexing::ExitCodes::UNEXPECTED_RESULT) { return UNEXPECTED_RESULT; }
        else { return UNKNOWN_ERROR; }
      }
    }

    e = switchScoreType_(peptide_ids, id_file_abs_path);
    if (e != EXECUTION_OK) return e;

    // TODO we could think about removing this limitation but it gets complicated quickly
    IDFilter::keepBestPeptideHits(peptide_ids, false); // strict = false

    // add to the (global) set of fixed and variable modifications
    const vector<std::string>& var_mods = protein_ids[0].getSearchParameters().variable_modifications;
    const vector<std::string>& fixed_mods = protein_ids[0].getSearchParameters().fixed_modifications;
    std::copy(var_mods.begin(), var_mods.end(), std::inserter(variable_modifications, variable_modifications.begin()));
    std::copy(fixed_mods.begin(), fixed_mods.end(), std::inserter(fixed_modifications, fixed_modifications.end()));

    // 'protein_references' is annotated by the PeptideIndexing run above (or already present in the
    // input) and is the only source of the theoretical uniqueness information that
    // IDFilter::keepUniquePeptidesPerProtein() needs in inferProteinGroups_(). Stripping it here
    // silently removed every peptide hit for '-protein_quantification strictly_unique_peptides'.
    const bool keep_protein_references = getStringOption_("protein_quantification") == "strictly_unique_peptides";

    // delete meta info to free some space
    for (PeptideIdentification& pid : peptide_ids)
    {
      // we currently can't clear the PeptideIdentification meta data
      // because the spectrum_reference is stored in the meta value (which it probably shouldn't)
      // TODO: pid.clearMetaInfo(); if we move it to the PeptideIdentification structure
      for (PeptideHit& ph : pid.getHits())
      {
        // TODO: we only have super inefficient meta value removal
        vector<std::string> keys;
        ph.getKeys(keys);
        for (const auto& k : keys)
        {
          if (!(StringUtils::hasSubstring(k, "_score")
            || StringUtils::hasSubstring(k, "q-value")
            || StringUtils::hasPrefix(k, "Luciphor_global_flr")
            || k == "target_decoy" // keep target_decoy information for QC
            || (keep_protein_references && k == "protein_references"))
            )
          {
            ph.removeMetaValue(k);
          }
        }
        // we only clear selected metavalues
        // ph.clearMetaInfo();
      }
    }

    ///////////////////////////////////////////////////////
    // annotate experimental design
    // check and reannotate mzML file in ID
    StringList id_msfile_ref;
    protein_ids[0].getPrimaryMSRunPath(id_msfile_ref);

    // fix other problems like missing MS run path annotations
    if (id_msfile_ref.empty())
    {
      OPENMS_LOG_WARN << "MS run path not set in ID file: " << id_file_abs_path << endl
                      << "Resetting reference to MS file provided at same input position." << endl;
    }
    else if (id_msfile_ref.size() == 1)
    {
      // Check if the annotated primary MS run filename matches the mzML filename (comparison by base name)
      const std::string& in_bn = File::stemName(mz_file_abs_path);
      const std::string& id_primaryMSRun_bn = File::stemName(id_msfile_ref[0]);

      if (in_bn != id_primaryMSRun_bn) // mismatch between annotation in ID file and provided mzML file
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
    id_msfile_ref = StringList {mz_file};
    protein_ids[0].setPrimaryMSRunPath(id_msfile_ref);
    protein_ids[0].setMetaValue("fraction_group", fraction_group);
    protein_ids[0].setMetaValue("fraction", fraction);

    // update identifiers to make them unique
    // fixes some bugs related to users splitting the original mzML and id files before running the analysis
    // in that case these files might have the same identifier
    const std::string old_identifier = protein_ids[0].getIdentifier();
    const std::string new_identifier = old_identifier + "_" + StringUtils::toStr(fraction_group) + "F" + StringUtils::toStr(fraction);
    protein_ids[0].setIdentifier(new_identifier);
    for (PeptideIdentification& p : peptide_ids)
    {
      if (p.getIdentifier() == old_identifier) { p.setIdentifier(new_identifier); }
      else { OPENMS_LOG_WARN << "Peptide ID identifier found not present in the protein ID" << endl; }
    }

    bool missing_spec_ref(false);
    for (const PeptideIdentification& pid : peptide_ids)
    {
      if (! pid.metaValueExists(Constants::UserParam::SPECTRUM_REFERENCE)
          || pid.getMetaValue(Constants::UserParam::SPECTRUM_REFERENCE).toString().empty())
      {
        missing_spec_ref = true;
        break;
      }
    }
    // reannotate spectrum references if missing
    if (missing_spec_ref)
    {
      OPENMS_LOG_WARN << "Warning: Identification file " << id_file_abs_path
                      << " contains IDs without meta value for the spectrum native id.\n"
                         "OpenMS will try to reannotate them by matching retention times between ID and spectra."
                      << endl;

      SpectrumMetaDataLookup::addMissingSpectrumReferences(peptide_ids, mz_file_abs_path, true);
    }

    // Deliberately last, i.e. after the spectrum-reference repair above: an identification with
    // no reference cannot be keyed at all, so reducing any earlier would silently skip exactly
    // the files whose references had to be reannotated.
    reduceToOnePerSpectrum_(peptide_ids, id_file_abs_path);

    return EXECUTION_OK;
  }

  /**
    @brief Reduce identifications this tool cannot tell apart to one per spectrum.

    ProteomicsLFQ measures one value per (spectrum, peptidoform, charge). Input that carries
    several identifications of that same triple - two search engines concatenated rather than
    combined, for instance - leaves it undecided which of them owns the measurement, and keeping
    all of them counts one identification several times in the PSM-level FDR and repeats it in
    every output. Keep the best-scoring one and say what went.
  **/
  void reduceToOnePerSpectrum_(PeptideIdentificationList& peptide_ids,
                               const std::string& id_file_abs_path)
  {
    const auto report = IDConflictResolverAlgorithm::reduceToOnePerSpectrum(peptide_ids);

    if (report.removed > 0)
    {
      OPENMS_LOG_WARN << "Warning: " << report.removed << " identification(s) in " << id_file_abs_path
                      << " repeat a (spectrum, peptidoform, charge) already claimed by another"
                         " identification, e.g. " << report.example << ".\n"
                         "Kept the best-scoring one of each group; the others would have counted"
                         " the same measurement more than once. To control this upstream, combine"
                         " search engine results with ConsensusID (-algorithm best"
                         " -keep_old_scores) rather than concatenating them." << endl;
    }
    if (report.inconsistent_score_direction > 0)
    {
      OPENMS_LOG_WARN << "Warning: " << report.inconsistent_score_direction << " group(s) of"
                         " repeated identifications in " << id_file_abs_path << " disagree on"
                         " whether a higher score is better, so no best one could be chosen."
                         " They were left as they are and will be counted more than once." << endl;
    }
    if (report.without_spectrum_reference > 0)
    {
      OPENMS_LOG_INFO << report.without_spectrum_reference << " identification(s) in "
                      << id_file_abs_path << " carry no spectrum reference and were therefore not"
                         " checked for repeats." << endl;
    }
    if (report.multiply_identified_spectra > 0)
    {
      OPENMS_LOG_INFO << report.multiply_identified_spectra << " spectra in " << id_file_abs_path
                      << " carry more than one identification naming different peptidoforms"
                         " (a chimeric spectrum, or search engines that disagree). Each is"
                         " quantified separately." << endl;
    }
  }


  void printMetaValues(const FeatureMap& tmp)
  {
    if (tmp.empty())
    {
      OPENMS_LOG_WARN << "printMetaValues called with empty FeatureMap.\n";
      return;
    }

    // extract meta value keys from the first element (which might be a normal or OffsetPeptide -> extract only the common ones)
    std::vector<std::string> keys;
    tmp[0].getKeys(keys);
    if (auto it = std::find(keys.begin(), keys.end(), "OffsetPeptide"); it != keys.end()) { keys.erase(it); }

    OPENMS_LOG_INFO << "keys: ";
    for (const auto& k : keys)
    {
      OPENMS_LOG_INFO << k << " ";
    }
    OPENMS_LOG_INFO << "\n";

    for (const auto& f : tmp)
    {
      if (f.metaValueExists("OffsetPeptide")) continue;
      for (const auto& k : keys)
      {
        OPENMS_LOG_INFO << f.getMetaValue(k) << " ";
      }
      OPENMS_LOG_INFO << "\n";
    }
  }

  /// Bumped when the stored layout changes in a way an older reader would misread.
  static constexpr int CHECKPOINT_SCHEMA_VERSION = 1;

  /**
    @brief What a checkpoint has to agree with before it may stand in for detecting a run.

    Assembled once per invocation. Everything here is compared, and a difference in any of it
    except @p openms_revision is fatal: a checkpoint that disagrees describes a different
    computation, and reusing it produces a result that looks fine and is not.
  */
  struct CheckpointContext
  {
    std::string fingerprint;   ///< the detect-relevant parameters, one "key=value" per line
    std::string design_row;    ///< fraction group, fraction, label and sample of this run
    std::string fasta_stamp;
    std::string openms_version;
    std::string openms_revision;
    // Paths rather than ready-made stamps: the stamp is taken at the point of comparison, so a
    // checkpoint rejected on cheaper grounds never causes the input to be looked at at all.
    // Empty when the file is not available at all, which is the case for a combining run.
    std::string mzml_path;
    std::string id_path;
  };

  /// FNV-1a. Not a security hash -- it identifies content, and unlike std::hash it is stable
  /// across builds and platforms, which a value written into a file has to be.
  static std::string digestString_(const std::string& s)
  {
    uint64_t h = 1469598103934665603ULL;
    for (unsigned char c : s) { h ^= c; h *= 1099511628211ULL; }
    std::ostringstream os;
    os << std::hex << std::setw(16) << std::setfill('0') << h;
    return os.str();
  }

  /**
    @brief Identify an input by its size and modification time, the way build tools do.

    Deliberately not a content hash. OpenMS's SHA-1 runs at about 135 MB/s -- measured, and some
    fifty times slower than simply reading the bytes -- so confirming that a 3 GB mzML has not
    changed costs more than a stat() by four orders of magnitude, on a path whose entire purpose
    is to avoid redoing work. make, ninja and ccache decide the same question the same way.

    Applied to every input, the FASTA included, so there is one rule rather than two. Databases
    reach 10 GB in metaproteomics and six-frame work, and the FASTA is checked by every detect job,
    so hashing it would add over a minute to each -- often more than the detection itself.

    The hole is the one those tools accept: a change preserving both size and modification time
    goes unnoticed. For the FASTA there is a second, sharper limitation, because it is the one
    input compared ACROSS machines: two copies of the same database staged separately have the
    same size but unrelated timestamps, so a run whose nodes each stage their own copy -- Nextflow
    in copy mode, containers, cloud executors -- finds every checkpoint disagreeing about a
    database that is in fact identical. This is built for the shared-filesystem workflow, where the
    FASTA is one file at one path and the comparison is exact; elsewhere, stage the database once
    and point every job at that copy.

    Directories (Bruker '.d', '.idparquet') are stamped from their sorted entries, since a
    directory's own size and timestamp say nothing about its contents.
  */
  std::string inputStamp_(const std::string& path) const
  {
    if (path.empty() || ! File::exists(path)) return "";

    auto stamp_of = [](const std::string& f) {
      return StringUtils::toStr(File::fileSize(f)) + ":" + StringUtils::toStr(File::getModificationTime(f));
    };

    if (! File::isDirectory(path)) return stamp_of(path);

    // Walked recursively, and keyed on the path relative to the root: a Bruker '.d' keeps content
    // in subdirectories, and File::fileList is a flat directory_iterator, so a flat listing would
    // be blind to most of one. Relative rather than basename, because two entries in different
    // subdirectories can share a name.
    std::vector<std::string> lines;
    std::error_code ec;
    const std::filesystem::path root(path);
    for (std::filesystem::recursive_directory_iterator it(root, ec), end; it != end && ! ec; it.increment(ec))
    {
      if (! it->is_regular_file(ec)) continue;
      lines.push_back(std::filesystem::relative(it->path(), root, ec).generic_string()
                      + "=" + stamp_of(it->path().string()));
    }
    if (ec) return ""; // unreadable: treat as unidentifiable rather than as "matches"
    std::sort(lines.begin(), lines.end());
    return "dir:" + digestString_(ListUtils::concatenate(lines, ";"));
  }

  /**
    @brief The parameters that shaped a run, as sorted "key=value" lines.

    Scoped by EXCLUDING what only the combining half uses, not by listing what the detecting half
    uses. The asymmetry is deliberate: forgetting to exclude a new combine-only parameter produces
    a loud false mismatch, while forgetting to include a new detection parameter in an allowlist
    would silently combine checkpoints that disagree.
  */
  std::string checkpointFingerprint_() const
  {
    static const std::vector<std::string> excluded_prefixes = {
      "Alignment:", "Linking:", "PipEcho:", "Protein Inference:", "ProteinQuantification:",
      "Posterior Error Probability:"};
    static const std::set<std::string> excluded_keys = {
      // Paths, not settings. 'fasta' in particular: the database is identified by its own stamp,
      // and leaving the path here would make two machines that mount the same database at
      // different points disagree about the *parameters* -- which is the normal case for the
      // distributed workflow this exists for, and would reject every checkpoint.
      "in", "ids", "design", "fasta", "out", "out_cxml", "out_msstats", "out_qpx",
      "feat_dir", "detect_only", "force_recompute",
      "threads", "debug", "log", "no_progress", "test", "force", "version", "write_ini", "ini",
      "proteinFDR", "psmFDR", "picked_proteinFDR", "FDR_type", "protein_inference",
      "protein_quantification", "alignment_order"};
    // Every entry above names a parameter this tool actually registers, except the command-line-only
    // ones ('version', 'write_ini', 'ini') and the two prefixes not yet registered as subsections
    // ('Protein Inference:', 'Posterior Error Probability:', both read by the combining half only).
    // Keep it that way: an entry matching nothing looks like it excludes something and does not, so
    // correcting its spelling later would drop a live setting out of the fingerprint without a word.

    std::vector<std::string> lines;
    for (auto it = getParam_().begin(); it != getParam_().end(); ++it)
    {
      const std::string key = it.getName();
      if (excluded_keys.count(key) != 0) continue;
      bool excluded = false;
      for (const auto& p : excluded_prefixes) { if (StringUtils::hasPrefix(key, p)) { excluded = true; break; } }
      if (excluded) continue;
      lines.push_back(key + "=" + it->value.toString());
    }
    std::sort(lines.begin(), lines.end());
    return ListUtils::concatenate(lines, "\n");
  }

  /// Name a checkpoint after the design row it belongs to, not after the file alone: two rows can
  /// share a basename in different fractions. The stem stays in front so a directory listing is
  /// still readable.
  static std::string checkpointName_(const std::string& mz_file, Size fraction_group, Size fraction)
  {
    return File::stemName(File::basename(mz_file)) + "_fg" + StringUtils::toStr(fraction_group)
           + "_f" + StringUtils::toStr(fraction) + ".featureParquet";
  }

  static std::string designRowKey_(Size fraction_group, Size fraction, const std::string& sample_name)
  {
    return "fg=" + StringUtils::toStr(fraction_group) + ";f=" + StringUtils::toStr(fraction)
           + ";sample=" + sample_name;
  }

  /**
    @brief Everything one MS run contributes to its fraction.

    This is the complete per-run result. Anything a run produces that is not in here leaves
    detectRun_ as a side effect instead, and side effects are precisely what does not survive
    a run being computed somewhere else. So this is a contract, not a convenience bundle: it
    is what a caller must be handed if it is to reproduce the run without re-reading the raw
    data.
  */
  struct RunDetection
  {
    FeatureMap features;                  ///< the quantified features of this run
    double fwhm = 0.0;                    ///< this run's chromatographic FWHM
    std::string id_ms_run_ref;            ///< primary MS run path as annotated in the ID file
    set<std::string> fixed_modifications; ///< search modifications seen in this run
    set<std::string> variable_modifications;
    std::string decoy_string;             ///< decoy affix inferred by PeptideIndexing
    bool decoy_prefix = true;
    bool decoy_inferred = false;          ///< false if no FASTA was given, so nothing was inferred
  };

  /// Report which keys of two fingerprints differ, so a mismatch names the setting rather than
  /// leaving an operator to diff INIs across nodes by hand.
  static std::string fingerprintDiff_(const std::string& mine, const std::string& theirs)
  {
    auto to_map = [](const std::string& s) {
      std::map<std::string, std::string> m;
      for (const auto& line : ListUtils::create<std::string>(s, '\n'))
      {
        const auto pos = line.find('=');
        if (pos != std::string::npos) m[line.substr(0, pos)] = line.substr(pos + 1);
      }
      return m;
    };
    const auto a = to_map(mine), b = to_map(theirs);
    std::vector<std::string> diffs;
    for (const auto& [k, v] : a)
    {
      auto it = b.find(k);
      if (it == b.end()) diffs.push_back("  " + k + ": '" + v + "' vs <absent in checkpoint>");
      else if (it->second != v) diffs.push_back("  " + k + ": '" + v + "' vs '" + it->second + "'");
    }
    for (const auto& [k, v] : b) { if (a.find(k) == a.end()) diffs.push_back("  " + k + ": <absent here> vs '" + v + "'"); }
    if (diffs.size() > 12) { diffs.resize(12); diffs.push_back("  ... and more"); }
    return ListUtils::concatenate(diffs, "\n");
  }

  /**
    @brief Write one run's checkpoint.

    Written under a temporary name inside @p dir and renamed into place only once every file is
    closed, so a checkpoint is either complete or absent and a resuming run can trust a plain
    existence check. The stamped meta values are removed again afterwards: the map continues into
    alignment and linking, and leaving them on it would put them in the consensusXML of a
    checkpointed run but not of a single-shot one.
  */
  bool writeCheckpoint_(const std::string& dir, const std::string& name,
                        RunDetection& rd, const CheckpointContext& ctx, bool replace_existing) const
  {
    FeatureMap& fm = rd.features;
    const std::vector<std::pair<std::string, DataValue>> stamped = {
      {"PLFQ:schema", DataValue(CHECKPOINT_SCHEMA_VERSION)},
      {"PLFQ:run_identifier", DataValue(rd.id_ms_run_ref.empty() ? "" : rd.id_ms_run_ref)},
      {"PLFQ:canonical_identifier", DataValue(fm.getProteinIdentifications().empty()
                                              ? std::string() : fm.getProteinIdentifications()[0].getIdentifier())},
      {"PLFQ:fwhm", DataValue(rd.fwhm)},
      {"PLFQ:decoy_string", DataValue(rd.decoy_string)},
      {"PLFQ:decoy_prefix", DataValue(rd.decoy_prefix ? 1 : 0)},
      {"PLFQ:decoy_inferred", DataValue(rd.decoy_inferred ? 1 : 0)},
      {"PLFQ:fixed_mods", DataValue(StringList(rd.fixed_modifications.begin(), rd.fixed_modifications.end()))},
      {"PLFQ:var_mods", DataValue(StringList(rd.variable_modifications.begin(), rd.variable_modifications.end()))},
      {"PLFQ:fingerprint", DataValue(ctx.fingerprint)},
      {"PLFQ:design_row", DataValue(ctx.design_row)},
      {"PLFQ:fasta_stamp", DataValue(ctx.fasta_stamp)},
      {"PLFQ:mzml_stamp", DataValue(inputStamp_(ctx.mzml_path))},
      {"PLFQ:id_stamp", DataValue(inputStamp_(ctx.id_path))},
      {"PLFQ:openms_version", DataValue(ctx.openms_version)},
      {"PLFQ:openms_revision", DataValue(ctx.openms_revision)}};
    for (const auto& [k, v] : stamped) { fm.setMetaValue(k, v); }

    // getUniqueName() rather than the pid: '-feat_dir' is documented as needing no coordination,
    // so two nodes may legitimately be detecting the same run into one shared directory, and a pid
    // is not unique between machines -- containerised executors hand out the same low pids by the
    // dozen. getUniqueName() folds in the hostname and a process-local counter as well.
    const std::string tmp = dir + "/." + name + ".tmp." + File::getUniqueName();
    if (File::exists(tmp)) { File::removeDirRecursively(tmp); }
    const bool ok = FeatureMapArrowIO::exportToParquet(fm, tmp);

    for (const auto& [k, v] : stamped) { fm.removeMetaValue(k); }

    if (! ok)
    {
      OPENMS_LOG_ERROR << "Could not write feature checkpoint for " << name << "\n";
      File::removeDirRecursively(tmp);
      return false;
    }

    const std::string final_path = dir + "/" + name;
    if (File::exists(final_path))
    {
      if (! replace_existing)
      {
        // Another process finished this run while we were writing it. Its checkpoint satisfies the
        // same contract, so keep it rather than replacing a file something may already be reading.
        OPENMS_LOG_INFO << "Checkpoint " << name << " appeared while it was being written; keeping the existing one.\n";
        File::removeDirRecursively(tmp);
        return true;
      }
      // -force_recompute: the whole point is to supersede what is there. Move the old one aside
      // first and delete it only once the new one is in place, so the name is never absent and a
      // failed rename cannot leave the run without a checkpoint at all.
      const std::string superseded = dir + "/." + name + ".superseded." + File::getUniqueName();
      File::removeDirRecursively(superseded);
      if (! File::rename(final_path, superseded, false))
      {
        OPENMS_LOG_ERROR << "Could not move the existing checkpoint aside: " << final_path << "\n";
        File::removeDirRecursively(tmp);
        return false;
      }
      if (! File::rename(tmp, final_path, false))
      {
        File::rename(superseded, final_path, false); // put the old one back
        File::removeDirRecursively(tmp);
        OPENMS_LOG_ERROR << "Could not move the recomputed checkpoint into place: " << final_path << "\n";
        return false;
      }
      File::removeDirRecursively(superseded);
      OPENMS_LOG_INFO << "Replaced checkpoint " << name << " (-force_recompute).\n";
      return true;
    }
    if (! File::rename(tmp, final_path, false))
    {
      // The existence check above is not atomic with the rename, so another invocation given the
      // same run can land its checkpoint in between. Detect-only invocations are documented as
      // needing no coordination, so losing that race has to be success, not an aborted job.
      File::removeDirRecursively(tmp);
      if (File::exists(final_path))
      {
        OPENMS_LOG_INFO << "Checkpoint " << name << " was written concurrently by another process; "
                        << "keeping theirs.\n";
        return true;
      }
      OPENMS_LOG_ERROR << "Could not move the feature checkpoint into place: " << final_path << "\n";
      return false;
    }
    return true;
  }

  enum class CheckpointState { ABSENT, USABLE, REJECTED };

  /// Load a checkpoint and decide whether it may stand in for detecting the run.
  CheckpointState loadCheckpoint_(const std::string& path, const CheckpointContext& ctx,
                                  RunDetection& rd, std::string& reason) const
  {
    if (! File::exists(path) || ! File::isDirectory(path)) return CheckpointState::ABSENT;

    for (const char* f : {"features.parquet", "psms.parquet", "proteins.parquet",
                          "protein_groups.parquet", "search_params.parquet"})
    {
      if (! File::exists(path + "/" + f))
      {
        reason = std::string("it is missing ") + f + ", so it is not a complete checkpoint";
        return CheckpointState::REJECTED;
      }
    }

    FeatureMap fm;
    if (! FeatureMapArrowIO::importFromParquet(path, fm))
    {
      reason = "it could not be read (see the error above)";
      return CheckpointState::REJECTED;
    }
    if (! fm.metaValueExists("PLFQ:schema"))
    {
      reason = "it carries no ProteomicsLFQ checkpoint metadata; it was not written by this tool";
      return CheckpointState::REJECTED;
    }
    if (static_cast<int>(fm.getMetaValue("PLFQ:schema")) != CHECKPOINT_SCHEMA_VERSION)
    {
      reason = "its checkpoint schema is version " + fm.getMetaValue("PLFQ:schema").toString()
               + ", this build writes version " + StringUtils::toStr(CHECKPOINT_SCHEMA_VERSION);
      return CheckpointState::REJECTED;
    }

    const auto stored = [&fm](const char* k) { return fm.metaValueExists(k) ? fm.getMetaValue(k).toString() : std::string(); };

    // Everything below reads these unconditionally, and getMetaValue throws on a missing key. A
    // schema-1 checkpoint always has them, so one that does not is malformed rather than merely
    // out of date -- report it as such instead of dying with an exception from the kernel classes.
    for (const char* k : {"PLFQ:fwhm", "PLFQ:decoy_prefix", "PLFQ:decoy_inferred",
                          "PLFQ:fixed_mods", "PLFQ:var_mods", "PLFQ:canonical_identifier"})
    {
      if (! fm.metaValueExists(k))
      {
        reason = std::string("it is missing the '") + k + "' entry, so it is malformed";
        return CheckpointState::REJECTED;
      }
    }

    if (stored("PLFQ:design_row") != ctx.design_row)
    {
      reason = "it was written for design row '" + stored("PLFQ:design_row") + "', not '" + ctx.design_row + "'";
      return CheckpointState::REJECTED;
    }
    if (stored("PLFQ:fingerprint") != ctx.fingerprint)
    {
      reason = "it was produced with different settings:\n" + fingerprintDiff_(ctx.fingerprint, stored("PLFQ:fingerprint"));
      return CheckpointState::REJECTED;
    }
    if (stored("PLFQ:fasta_stamp") != ctx.fasta_stamp)
    {
      reason = "the FASTA does not match the one it was indexed against. These are compared by size "
               "and modification time, so this also fires when the same database was staged separately "
               "for the run that wrote the checkpoint -- point every job at one copy on shared storage, "
               "or re-detect with -force_recompute";
      return CheckpointState::REJECTED;
    }
    if (stored("PLFQ:openms_version") != ctx.openms_version
        || stored("PLFQ:openms_revision") != ctx.openms_revision)
    {
      // No override. Detection behaviour can change between any two builds, and nothing here can
      // tell whether it did, so combining checkpoints across builds would be a guess dressed up as
      // a result. The intended use -- one build, many machines -- never reaches this, and
      // -force_recompute is always available when it does.
      reason = "it was written by OpenMS " + stored("PLFQ:openms_version") + " ("
               + stored("PLFQ:openms_revision") + "), this is " + ctx.openms_version + " ("
               + ctx.openms_revision + "). Detect the run again with -force_recompute, or delete "
               "the checkpoint directory";
      return CheckpointState::REJECTED;
    }

    // Last, because these are the only gates that touch the filesystem rather than the metadata
    // already in hand -- a stat() per input, and for a directory input one per entry in it. A
    // checkpoint that disagrees about anything cheaper is rejected before reaching them.
    // Only checkable when the inputs are at hand; a combining run has neither.
    if (! ctx.mzml_path.empty() && stored("PLFQ:mzml_stamp") != inputStamp_(ctx.mzml_path))
    {
      reason = "the spectra file has changed since it was written";
      return CheckpointState::REJECTED;
    }
    if (! ctx.id_path.empty() && stored("PLFQ:id_stamp") != inputStamp_(ctx.id_path))
    {
      reason = "the identification file has changed since it was written";
      return CheckpointState::REJECTED;
    }

    rd.fwhm = fm.getMetaValue("PLFQ:fwhm");
    rd.decoy_string = stored("PLFQ:decoy_string");
    rd.decoy_prefix = static_cast<int>(fm.getMetaValue("PLFQ:decoy_prefix")) != 0;
    rd.decoy_inferred = static_cast<int>(fm.getMetaValue("PLFQ:decoy_inferred")) != 0;
    rd.id_ms_run_ref = stored("PLFQ:run_identifier");
    const StringList fixed = fm.getMetaValue("PLFQ:fixed_mods");
    const StringList var = fm.getMetaValue("PLFQ:var_mods");
    rd.fixed_modifications.insert(fixed.begin(), fixed.end());
    rd.variable_modifications.insert(var.begin(), var.end());

    // The Parquet reader synthesizes a fresh run identifier on load and treats the stored one as
    // informational, so the canonical identifier has to be put back by hand -- on the protein run
    // and on every peptide identification that references it.
    const std::string canonical = stored("PLFQ:canonical_identifier");
    if (! canonical.empty() && ! fm.getProteinIdentifications().empty())
    {
      fm.getProteinIdentifications()[0].setIdentifier(canonical);
      for (auto& f : fm) { for (auto& p : f.getPeptideIdentifications()) { p.setIdentifier(canonical); } }
      for (auto& p : fm.getUnassignedPeptideIdentifications()) { p.setIdentifier(canonical); }
    }

    for (const auto& k : {"PLFQ:schema", "PLFQ:run_identifier", "PLFQ:canonical_identifier", "PLFQ:fwhm",
                          "PLFQ:decoy_string", "PLFQ:decoy_prefix", "PLFQ:decoy_inferred", "PLFQ:fixed_mods",
                          "PLFQ:var_mods", "PLFQ:fingerprint", "PLFQ:design_row", "PLFQ:fasta_stamp",
                          "PLFQ:mzml_stamp", "PLFQ:id_stamp", "PLFQ:openms_version", "PLFQ:openms_revision"})
    {
      fm.removeMetaValue(k);
    }

    // The Parquet reader annotates every PeptideHit with 'scan' and 'reference_file_name', derived
    // from columns it stores for analytics, whether or not the written hit carried them. A detected
    // run never carries them -- loadAndCleanupIDFile_ strips every PeptideHit meta value that is not
    // a score, a q-value, a Luciphor FLR, target_decoy or (for strictly unique peptides)
    // protein_references -- so leaving them on would make a checkpointed run emit two mzTab opt_
    // columns a single-shot run does not, and the two would stop being comparable. Restore what was
    // written.
    auto drop_invented = [](PeptideIdentificationList& pids) {
      for (auto& pid : pids)
      {
        for (auto& hit : pid.getHits())
        {
          hit.removeMetaValue("scan");
          hit.removeMetaValue("reference_file_name");
        }
      }
    };
    for (auto& f : fm) { drop_invented(f.getPeptideIdentifications()); }
    drop_invented(fm.getUnassignedPeptideIdentifications());

    rd.features = std::move(fm);
    return CheckpointState::USABLE;
  }

  /**
    @brief Detect and quantify the features of a single MS run.

    Depends on nothing but its own inputs, which is what allows the runs of a fraction to be
    processed in any order (and, eventually, anywhere). Alignment, linking, inference and
    quantification need every run at once and stay in quantifyFraction_ and main_.
  */
  ExitCodes detectRun_(
    const std::string& mz_file,
    const std::string& id_file_abs_path,
    const std::string& in_db,
    const Size fraction,
    const Size fraction_group,
    RunDetection& out)
  {
    writeDebug_("Processing file: " + mz_file, 1);

    vector<ProteinIdentification> protein_ids;
    PeptideIdentificationList peptide_ids;

    {
      ExitCodes e = loadAndCleanupIDFile_(id_file_abs_path, mz_file, in_db, fraction_group, fraction, protein_ids, peptide_ids, out.fixed_modifications,
                                          out.variable_modifications, out.decoy_string, out.decoy_prefix, out.decoy_inferred);
      if (e != EXECUTION_OK) return e;
    }

    MSExperiment ms_centroided;
    bool is_im_peak_data = false;
    const bool mass_recalibration = (getStringOption_("mass_recalibration") == "true");

    // Chromatographic FWHM of THIS run. Scoped to the iteration so that a run which cannot
    // measure one (IM_PEAK without seeds) falls back to the documented default instead of
    // silently inheriting the previous run's value.
    double median_fwhm = 0.0;

    {
      ExitCodes e = loadAndPreprocess_(mz_file, ms_centroided, is_im_peak_data);
      if (e != EXECUTION_OK) { return e; }
    }

    SpectrumMetaDataLookup::addMissingFAIMSToPeptideIDs(peptide_ids, ms_centroided);

    if (is_im_peak_data && mass_recalibration)
    {
      OPENMS_LOG_WARN << "Warning: mass_recalibration is not supported for .d input. Disabling.\n";
    }

    if (mass_recalibration && !is_im_peak_data)
    {
      std::string debug_output_basename = (debug_level_ > 666) ? id_file_abs_path : "";
      DDAWorkflowCommons::recalibrateMS1(ms_centroided, peptide_ids, debug_output_basename);
    }

    if (!is_im_peak_data)
    {
      median_fwhm = DDAWorkflowCommons::estimateMedianChromatographicFWHM(ms_centroided);
      OPENMS_LOG_INFO << "Median chromatographic FWHM: " << median_fwhm << "\n";
    }
    // For .d/IM_PEAK: FWHM is estimated from Biosaur2 features below (after seed generation)

    // For .d/IM_PEAK + targeted_only: no Biosaur2 seeds to estimate FWHM from, use default
    if (is_im_peak_data && median_fwhm == 0.0)
    {
      median_fwhm = 30.0;
      OPENMS_LOG_INFO << "Using default FWHM of " << median_fwhm << " s for .d input (no seed-based estimate available).\n";
    }

    StringList id_msfile_ref;
    protein_ids[0].getPrimaryMSRunPath(id_msfile_ref);
    out.id_ms_run_ref = id_msfile_ref[0];

    FeatureMap seeds;
    seeds.setPrimaryMSRunPath({mz_file});

    const bool targeted_only = getStringOption_("targeted_only") != "false";

    if (! targeted_only)
    {
      std::string seeding_algorithm = getStringOption_("Seeding:algorithm");

      // Force biosaur2 for IM_PEAK data (Multiplex cannot handle it)
      if (is_im_peak_data && seeding_algorithm != "biosaur2")
      {
        OPENMS_LOG_WARN << "Warning: IM_PEAK data detected. Forcing Seeding:algorithm to 'biosaur2' "
                        << "(FeatureFinderMultiplex does not support IM_PEAK format).\n";
        seeding_algorithm = "biosaur2";
      }

      if (seeding_algorithm == "biosaur2")
      {
        OPENMS_LOG_INFO << "Using Biosaur2Algorithm for seed detection.\n";
        Biosaur2Algorithm bio;
        Param bio_param = getParam_().copy("Seeding:Biosaur2:", true);
        bio.setParameters(bio_param);
        // Biosaur2 leaves the processed spectra in ms_data_ on both its paths
        // (in-place for non-FAIMS, reassembled after the per-CV split for FAIMS),
        // so we can move here and retrieve the data after run().
        bio.setMSData(std::move(ms_centroided));

        bio.run(seeds);
        OPENMS_LOG_INFO << "Biosaur2 produced " << seeds.size() << " seed features.\n";

        // For .d/IM_PEAK: estimate FWHM from Biosaur2 feature convex hulls
        // Use FWHM ≈ 0.59 * (rt_end - rt_start), the Gaussian correction from base width to half-max width
        if (is_im_peak_data)
        {
          std::vector<double> fwhm_values;
          fwhm_values.reserve(seeds.size());
          for (const auto& f : seeds)
          {
            const auto& hulls = f.getConvexHulls();
            if (hulls.empty()) continue;
            const auto& bb = hulls[0].getBoundingBox();
            double base_width = bb.maxX() - bb.minX(); // full RT extent
            double fwhm = base_width * 0.59;           // Gaussian: FWHM ≈ 2.35σ, base ≈ 4σ → ratio ≈ 0.59
            if (fwhm > 0.0) { fwhm_values.push_back(fwhm); }
          }
          if (fwhm_values.size() >= 10)
          {
            median_fwhm = Math::median(fwhm_values.begin(), fwhm_values.end());
          }
          else
          {
            median_fwhm = 30.0; // fallback
            OPENMS_LOG_WARN << "Warning: Too few Biosaur2 features (" << fwhm_values.size()
                            << ") to estimate FWHM. Using default: " << median_fwhm << " seconds.\n";
          }
          OPENMS_LOG_INFO << "Median chromatographic FWHM (from Biosaur2): " << median_fwhm << "\n";
        }

        // Retrieve ms_data_ back from Biosaur2 (preserved on both the FAIMS and non-FAIMS paths)
        ms_centroided = std::move(bio.getMSData());
      }
      else
      {
        DDAWorkflowCommons::calculateSeeds(ms_centroided, getDoubleOption_("Seeding:intThreshold"), seeds, median_fwhm, 2, 5);
      }

      if (debug_level_ > 666)
      {
        FileHandler().storeFeatures("debug_seeds_fraction_" + StringUtils::toStr(fraction) + "_" + StringUtils::toStr(fraction_group) + ".featureXML", seeds, {FileTypes::FEATUREXML}, log_type_);
      }
    }

    // median_fwhm is final for this run: measured from the raw data, from the Biosaur2
    // seeds, or defaulted. Record it for the fraction-level linking window below.
    out.fwhm = median_fwhm;

    /////////////////////////////////////////////////
    // Feature detection
    FeatureMap fm;

    // Run FeatureFinderIdentification
    FeatureFinderIdentificationAlgorithm ffi;
    ffi.getMSData().swap(ms_centroided);
    ffi.getProgressLogger().setLogType(log_type_);

    Param ffi_param = getParam_().copy("PeptideQuantification:", true);
    ffi_param.setValue("detect:peak_width", 5.0 * median_fwhm);
    ffi_param.setValue("debug", debug_level_); // pass down debug level

    // PIP-ECHO's isotope-distribution MBR feature correlates the acceptor's
    // observed isotope envelope against the donor peptide's theoretical
    // envelope; the default 2 isotopes give a degenerate Pearson. Extract 3
    // (enough for a meaningful correlation) only on the pip_echo path. More
    // (e.g. 5) measurably perturbs FFID's extraction/candidate set and costs
    // direct quantifications, so 3 is the sweet spot. The default linker path
    // (and its TOPP reference outputs) is untouched.
    if (getStringOption_("pip_echo") == "true")
    {
      ffi_param.setValue("extract:n_isotopes", 3);
    }

    // Note: no IM_window override needed — BrukerTimsFile loads with built-in IM centroiding
    // (ms1_centroid_mz_ppm=5, ms1_centroid_im_pct=3), so data is IM_CENTROIDED and the
    // default IM_window=0.06 is appropriate for matching centroided IM positions.

    double feature_with_id_min_score = getDoubleOption_("feature_with_id_min_score");
    double feature_without_id_min_score = getDoubleOption_("feature_without_id_min_score");
    const bool filter_by_quant_scores = (feature_with_id_min_score > 0.0) && (targeted_only || (feature_without_id_min_score > 0.0));

    if (filter_by_quant_scores)
    {
      OPENMS_LOG_INFO << "Adding offset peptides as quant. decoys.\n";
      ffi_param.setValue("add_mass_offset_peptides", 10.005);
    }

    ffi.setParameters(ffi_param);
    writeDebug_("Parameters passed to FeatureFinderIdentification algorithm", ffi_param, 3);

    {
      vector<ProteinIdentification> ext_protein_ids;
      PeptideIdentificationList ext_peptide_ids;

      ffi.run(peptide_ids, protein_ids, ext_peptide_ids, ext_protein_ids,
              fm, // fills fm
              seeds, mz_file);
    }

    if (filter_by_quant_scores)
    {
      SimpleSVM::PredictorMap predictors;
      map<Size, double> labels;
      size_t current_row = 0;
      size_t quant_target {}, quant_decoy {};

      Math::RandomShuffler shuffler;
      std::vector<size_t> randomized_indices(fm.size());
      std::iota(randomized_indices.begin(), randomized_indices.end(), 0);

      for (auto& i : randomized_indices)
      {
        const auto& f = fm[i];
        predictors["var_library_sangle"].push_back(f.getMetaValue("var_library_sangle"));
        predictors["var_xcorr_shape"].push_back(f.getMetaValue("var_xcorr_shape"));
        predictors["total_xic"].push_back(f.getMetaValue("total_xic"));
        predictors["var_elution_model_fit_score"].push_back(f.getMetaValue("var_elution_model_fit_score"));

        bool is_offset = f.metaValueExists("OffsetPeptide");
        bool has_id = ! f.getPeptideIdentifications().empty();
        if (is_offset)
        {
          if (quant_decoy < 1000)
          {
            labels[current_row] = 0.0;
            ++quant_decoy;
          }
        }
        else
        {
          if (has_id && quant_target < 1000)
          {
            labels[current_row] = 1.0;
            ++quant_target;
          }
        }
        ++current_row;
      }

      if (quant_decoy > 4 && quant_target > 4)
      {
        SimpleSVM svm;
        Param svm_param = svm.getParameters();
        svm_param.setValue("kernel", "linear");
        svm_param.setValue("log2_C", ListUtils::create<double>("-5,-1,1,5,7,11,15"));
        svm_param.setValue("log2_p", ListUtils::create<double>("-15,-9,-6,-3.32192809489,0,3.32192809489,6,9,15"));

        svm.setParameters(svm_param);
        svm.setup(predictors, labels);
        vector<SimpleSVM::Prediction> predictions;
        OPENMS_LOG_INFO << "Predicting class probabilities:" << endl;
        svm.predict(predictions);
        std::map<std::string, double> feature_weights;
        svm.getFeatureWeights(feature_weights);

        size_t current_row_pred {};
        for (auto& i : randomized_indices)
        {
          auto& f = fm[i];
          f.setMetaValue("p_quant", (double)predictions[current_row_pred].probabilities[1]);
          ++current_row_pred;
        }
      }

      if (quant_decoy > 4 && quant_target > 4)
      {
        fm.erase(std::remove_if(fm.begin(), fm.end(),
                                [&](const Feature& f) {
                                  double quant_score = f.getMetaValue("p_quant");
                                  bool is_offset = f.metaValueExists("OffsetPeptide");
                                  bool has_id = ! f.getPeptideIdentifications().empty();
                                  bool untargeted_feature = ! is_offset && ! has_id;
                                  bool is_feature_with_id = ! is_offset && has_id;

                                  if (is_feature_with_id && quant_score < feature_with_id_min_score) return true;
                                  if (untargeted_feature && quant_score < feature_without_id_min_score) return true;
                                  if (is_offset && quant_score < feature_without_id_min_score) return true;
                                  return false;
                                }),
                 fm.end());

        fm.erase(std::remove_if(fm.begin(), fm.end(), [](const Feature& f) { return f.metaValueExists("OffsetPeptide"); }), fm.end());
      }
    }

    // "IM" and "n_scans" are written by Biosaur2 (ion-mobility value and true
    // MS1 scan count); kept for forward compatibility with a direct-Biosaur2
    // feature path (PIP-ECHO ion-mobility scoring reads "IM" as a fallback to
    // "IM_median"; "n_scans" feeds the planned scan-count feature, see #9655).
    // FFID output features carry neither key, so this is a no-op for the
    // current FFID-based path.
    // "IM"/"n_scans" carry no payload on the FFID path, so keeping them is a no-op.
    unordered_set<std::string> keep_meta = {"OffsetPeptide", "IM_median", "IM_min", "IM_max", "IM", "n_scans"};
    // "masserror_ppm" is FFID's observed mass-trace deviation from the peptide's
    // theoretical m/z (DIA mass-difference score) -- the only real per-feature
    // mass-accuracy signal (FFID seeds the feature's *position* m/z to theoretical,
    // so PIP-ECHO cannot recover the error from getMZ()). "pipecho_obs_envelope" is
    // the observed isotope envelope (subordinate intensities), snapshotted below
    // before subordinates are cleared, for the isotope-distribution feature. BOTH
    // are needed ONLY by PIP-ECHO; keeping/writing them on the default QT-linker path
    // leaks them into the consensusXML/mzTab output (shifting TOPP_ProteomicsLFQ
    // references), so gate on pip_echo.
    const bool pip_echo_enabled = getStringOption_("pip_echo") == "true";
    if (pip_echo_enabled) { keep_meta.insert("masserror_ppm"); keep_meta.insert("pipecho_obs_envelope"); }
    for (auto & f : fm)
    {
      std::vector<std::string> keys;
      f.getKeys(keys);
      for (const auto& k : keys)
      {
        if (auto it = keep_meta.find(k); it == keep_meta.end()) f.removeMetaValue(k);
      }
      // Snapshot the observed isotope envelope before subordinates are dropped
      // (pip_echo only -- otherwise it would leak into the default-path output,
      // since this writes AFTER the meta-strip above).
      if (pip_echo_enabled && ! f.getSubordinates().empty())
      {
        DoubleList envelope;
        envelope.reserve(f.getSubordinates().size());
        for (const Feature& sub : f.getSubordinates()) { envelope.push_back(sub.getIntensity()); }
        f.setMetaValue("pipecho_obs_envelope", envelope);
      }
      f.setSubordinates({});
      f.setConvexHulls({});
    }

    IDConflictResolverAlgorithm::resolve(fm, getStringOption_("keep_feature_top_psm_only") == "false");

    out.features = std::move(fm);

    if (debug_level_ > 666)
    {
      FileHandler().storeFeatures("debug_fraction_" + StringUtils::toStr(fraction) + "_" + StringUtils::toStr(fraction_group) + ".featureXML", out.features, {FileTypes::FEATUREXML}, log_type_);
    }

    return EXECUTION_OK;
  }

  ExitCodes quantifyFraction_(
    const pair<UInt, std::vector<std::string> > & ms_files,
    const map<std::string, std::string>& mzfile2idfile,
    const std::string& in_db,
    ConsensusMap & consensus_fraction,
    set<std::string>& fixed_modifications,
    set<std::string>& variable_modifications)
  {
    vector<TransformationDescription> transformations;
    vector<FeatureMap> feature_maps;
    const Size fraction = ms_files.first;

    writeDebug_("Processing fraction number: " + StringUtils::toStr(fraction) + "\nFiles: ",  1);
    for (std::string const & mz_file : ms_files.second) { writeDebug_(mz_file,  1); }

    StringList id_MS_run_ref;
    StringList in_MS_run = ms_files.second;
    const auto& path_label_to_fractiongroup = design_.getPathLabelToFractionGroupMapping(true);
    // By value, and hoisted: both accessors build and return a fresh map, so calling one per run
    // would rebuild it for every run, and binding a reference to an element of that temporary
    // would leave it dangling at the end of the statement.
    const auto path_label_to_sample = design_.getPathLabelToSampleMapping(true);

    // One chromatographic FWHM per MS run of this fraction. The linker's RT tolerance is a
    // property of the fraction, so it is derived from all of them (below) rather than from
    // whichever run the loop happened to end on.
    std::vector<double> run_fwhms;
    run_fwhms.reserve(ms_files.second.size());

    for (std::string const & mz_file : ms_files.second)
    {
      const Size fraction_group = path_label_to_fractiongroup.at({File::basename(mz_file), 1});
      const unsigned sample_idx = path_label_to_sample.at({File::basename(mz_file), 1});
      const std::string design_row = designRowKey_(fraction_group, fraction,
                                                   design_.getSampleSection().getSampleName(sample_idx));

      // Is this run's identification file available? A combining invocation has neither -in nor
      // -ids and relies entirely on checkpoints, so this is allowed to be absent.
      const auto id_it = mzfile2idfile.find(File::absolutePath(mz_file));
      const std::string id_file_abs_path =
        (id_it == mzfile2idfile.end()) ? std::string() : File::absolutePath(id_it->second);

      RunDetection rd;
      bool from_checkpoint = false;

      if (! feat_dir_.empty())
      {
        CheckpointContext ctx = checkpoint_ctx_;
        ctx.design_row = design_row;
        // Only stampable when the inputs are at hand. A combining run cannot notice that an mzML
        // changed after its checkpoint was written; that is the resuming run's job.
        if (! id_file_abs_path.empty())
        {
          ctx.mzml_path = File::absolutePath(mz_file);
          ctx.id_path = id_file_abs_path;
        }

        const std::string ckpt = feat_dir_ + "/" + checkpointName_(mz_file, fraction_group, fraction);
        if (getFlag_("force_recompute"))
        {
          if (File::exists(ckpt)) { OPENMS_LOG_INFO << "Recomputing " << File::basename(mz_file) << " (-force_recompute).\n"; }
        }
        else
        {
          std::string why;
          switch (loadCheckpoint_(ckpt, ctx, rd, why))
          {
            case CheckpointState::USABLE:
              OPENMS_LOG_INFO << "Reusing feature checkpoint for " << File::basename(mz_file) << ".\n";
              from_checkpoint = true;
              ++checkpoints_reused_;
              break;
            case CheckpointState::REJECTED:
              // Never silently recompute over a checkpoint that disagrees: the disagreement is the
              // information. Recomputing would hide that the operator combined an inconsistent set.
              OPENMS_LOG_FATAL_ERROR << "Feature checkpoint '" << ckpt << "' cannot be used because "
                                     << why << ".\nDelete it, or pass -force_recompute to detect this "
                                     << "run again and overwrite it.\n";
              return INCOMPATIBLE_INPUT_DATA;
            case CheckpointState::ABSENT:
              break;
          }
        }
      }

      if (! from_checkpoint)
      {
        if (id_file_abs_path.empty())
        {
          OPENMS_LOG_FATAL_ERROR << "No feature checkpoint for '" << File::basename(mz_file)
                                 << "' and no identification file given for it. A combining run needs a "
                                 << "checkpoint for every run of the design; otherwise pass -in and -ids "
                                 << "so it can be detected.\n";
          return INCOMPATIBLE_INPUT_DATA;
        }
        // Stamp the inputs before reading them, not after. Detection takes minutes; stamping
        // afterwards would pair features derived from the old content with a stamp describing the
        // new, and a later resume would accept that pairing as current.
        const std::string mzml_before = feat_dir_.empty() ? "" : inputStamp_(File::absolutePath(mz_file));
        const std::string id_before = feat_dir_.empty() ? "" : inputStamp_(id_file_abs_path);

        ExitCodes e = detectRun_(mz_file, id_file_abs_path, in_db, fraction, fraction_group, rd);
        if (e != EXECUTION_OK) { return e; }
        ++runs_detected_;

        if (! feat_dir_.empty())
        {
          if (inputStamp_(File::absolutePath(mz_file)) != mzml_before
              || inputStamp_(id_file_abs_path) != id_before)
          {
            OPENMS_LOG_FATAL_ERROR << "The inputs for '" << File::basename(mz_file)
                                   << "' changed while it was being processed, so no checkpoint can "
                                   << "honestly describe this run. Nothing was written for it.\n";
            return INCOMPATIBLE_INPUT_DATA;
          }
          CheckpointContext ctx = checkpoint_ctx_;
          ctx.design_row = design_row;
          ctx.mzml_path = File::absolutePath(mz_file);
          ctx.id_path = id_file_abs_path;
          if (! writeCheckpoint_(feat_dir_, checkpointName_(mz_file, fraction_group, fraction), rd, ctx,
                                 getFlag_("force_recompute")))
          {
            return CANNOT_WRITE_OUTPUT_FILE;
          }
        }
      }

      // Fold the run's contribution into the fraction. Every one of these was previously an
      // in-place write from inside the loop body; making them explicit here is the point of
      // the split.
      fixed_modifications.insert(rd.fixed_modifications.begin(), rd.fixed_modifications.end());
      variable_modifications.insert(rd.variable_modifications.begin(), rd.variable_modifications.end());
      recordDecoyAffix_(rd.decoy_string, rd.decoy_prefix, rd.decoy_inferred,
                        id_file_abs_path.empty() ? ("checkpoint of " + File::basename(mz_file)) : id_file_abs_path);
      id_MS_run_ref.push_back(rd.id_ms_run_ref);
      run_fwhms.push_back(rd.fwhm);
      feature_maps.emplace_back(std::move(rd.features));
    }

    if (getFlag_("detect_only"))
    {
      // Everything below needs every run of the fraction at once. A detect-only invocation exists
      // to produce the checkpoints and stop, so that the runs can be spread over machines and
      // combined later by an invocation that has no raw data at all.
      return EXECUTION_OK;
    }

    auto validation_result = File::validateMatchingFileNames(in_MS_run, id_MS_run_ref, true, true);
    switch (validation_result)
    {
      case File::MatchingFileListsStatus::SET_MISMATCH:
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "MS run path reference in ID files and spectra filenames differ.");
      case File::MatchingFileListsStatus::ORDER_MISMATCH:
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "MS run path reference in ID files and spectra filenames match but order differs.");
      case File::MatchingFileListsStatus::MATCH:
        writeLogInfo_("ID files reference same names as spectra files.");
        break;
    }

    // The linker derives its RT tolerance from this (2*max_alignment_diff + 2*fraction_fwhm,
    // and PIP-ECHO's local_rt window). It describes the fraction as a whole, so take the median
    // over its runs. Previously this read whatever the per-run variable happened to hold after
    // the loop, i.e. the last run's value, which made the tolerance depend on input order.
    double fraction_fwhm = 0.0;
    if (!run_fwhms.empty())
    {
      fraction_fwhm = Math::median(run_fwhms.begin(), run_fwhms.end());
      if (run_fwhms.size() > 1)
      {
        const auto minmax = std::minmax_element(run_fwhms.begin(), run_fwhms.end());
        OPENMS_LOG_INFO << "Median chromatographic FWHM across the " << run_fwhms.size()
                        << " runs of this fraction: " << fraction_fwhm
                        << " s (range " << *minmax.first << " - " << *minmax.second << " s)\n";
      }
    }

    alignAndLink_(feature_maps, consensus_fraction, transformations, fraction_fwhm);

    if (feature_maps.size() > 1)
    {
      addDataProcessing_(consensus_fraction, getProcessingInfo_(DataProcessing::ALIGNMENT));
      addDataProcessing_(consensus_fraction, getProcessingInfo_(DataProcessing::FEATURE_GROUPING));
    }

    Size j(0);
    const auto& path_label_to_sampleidx = design_.getPathLabelToSampleMapping(true);
    for (std::string const & mz_file : ms_files.second)
    {
      const Size curr_fraction_group = path_label_to_fractiongroup.at({File::basename(mz_file), 1});
      consensus_fraction.getColumnHeaders()[j].label = "label-free";
      consensus_fraction.getColumnHeaders()[j].filename = mz_file;
      consensus_fraction.getColumnHeaders()[j].unique_id = feature_maps[j].getUniqueId();
      consensus_fraction.getColumnHeaders()[j].setMetaValue("fraction", fraction);
      consensus_fraction.getColumnHeaders()[j].setMetaValue("fraction_group", curr_fraction_group);
      const auto& sample_index = path_label_to_sampleidx.at({File::basename(mz_file), 1});
      const auto& sample_name = design_.getSampleSection().getSampleName(sample_index);
      consensus_fraction.getColumnHeaders()[j].setMetaValue("sample_name", sample_name);
      ++j;
    }

    consensus_fraction.applyMemberFunction(&UniqueIdInterface::setUniqueId);
    consensus_fraction.sortPeptideIdentificationsByMapIndex();
    IDConflictResolverAlgorithm::resolve(consensus_fraction, true);

    // No normalization here, deliberately. A median scaling used to run at this point unless
    // -out_msstats was requested, which meant the same input written to the same
    // -out_qpx directory held normalized or raw intensities depending on whether an unrelated
    // output file was also asked for - and neither the Parquet nor the mzTab recorded which.
    // Normalizing per fraction and then summing across fractions (PeptideAndProteinQuant.cpp) is
    // also not well defined: it is only ratio-preserving when the reference is the same sample in
    // every fraction, which a non-rectangular design does not guarantee. See the tool
    // documentation above for the three supported ways to normalize.

    return EXECUTION_OK;
  }


  ExitCodes inferProteinGroups_(ConsensusMap& consensus,
    const set<std::string>& fixed_modifications)
  {
    // Note: protein sequences (for coverage) and the 'protein_references' uniqueness annotation used below
    // come from the PeptideIndexing run in loadAndCleanupIDFile_(); there is no second indexing pass here.

    //-------------------------------------------------------------
    // Protein inference
    //-------------------------------------------------------------
    // TODO: This needs to be rewritten to work directly on the quant data.
    //  of course we need to provide options to keep decoys and unassigned PSMs all the way through quantification.
    // TODO: Think about ProteinInference on IDs only merged per condition
    bool groups = getStringOption_("protein_quantification") != "strictly_unique_peptides";
    bool bayesian = getStringOption_("protein_inference") == "bayesian";
    bool greedy_group_resolution = getStringOption_("protein_quantification") == "shared_peptides";

    // Study-wide inference operates on a single merged ID run.
    ConsensusMapMergerAlgorithm cmerge;
    // The following will result in a SINGLE protein run for the whole consensusMap.
    cmerge.mergeAllIDRuns(consensus);

    if (! bayesian) // simple aggregation
    {
      BasicProteinInferenceAlgorithm bpia;
      auto bpiaparams = bpia.getParameters();
      bpiaparams.setValue("annotate_indistinguishable_groups", groups ? "true" : "false");
      bpiaparams.setValue("greedy_group_resolution", greedy_group_resolution ? "true" : "false");
      bpia.setParameters(bpiaparams);

      // TODO parameterize if unassigned IDs without feature should contribute?
      bpia.run(consensus, consensus.getProteinIdentifications()[0], true);
    }
    else // if (bayesian)
    {
      BayesianProteinInferenceAlgorithm bayes;
      auto bayesparams = bayes.getParameters();
      // We need all PSMs to collect all possible modifications, to do spectral counting and to do PSM FDR.
      // In theory, if none is needed we can save memory. For quantification,
      // we basically discard peptide+PSM information from inference and use the info from the cMaps.
      bayesparams.setValue("keep_best_PSM_only", "false");
      bayes.setParameters(bayesparams);
      // bayesian inference automatically annotates groups, therefore remove them later
      bayes.inferPosteriorProbabilities(consensus, greedy_group_resolution);
      if (! groups)
      {
        // should be enough to just clear the groups. Only indistinguishable will be annotated above.
        consensus.getProteinIdentifications()[0].getIndistinguishableProteins().clear();
      }
    }

    // TODO think about order of greedy resolution, FDR calc and filtering

    //-------------------------------------------------------------
    // Protein (and additional peptide?) FDR
    //-------------------------------------------------------------
    const double max_fdr = getDoubleOption_("proteinFDR");
    const bool picked = getStringOption_("picked_proteinFDR") == "true";

    // TODO use new FDR_type parameter
    const double max_psm_fdr = getDoubleOption_("psmFDR");
    FalseDiscoveryRate fdr;
    if (getFlag_("PeptideQuantification:quantify_decoys"))
    {
      Param fdr_param = fdr.getParameters();
      fdr_param.setValue("add_decoy_peptides", "true");
      fdr_param.setValue("add_decoy_proteins", "true");
      fdr.setParameters(fdr_param);
    }

    // ensure that only one final inference result is generated for now
    assert(consensus.getProteinIdentifications().size() == 1);

    auto& overall_proteins = consensus.getProteinIdentifications()[0];
    if (! picked) { fdr.applyBasic(overall_proteins); }
    else { fdr.applyPickedProteinFDR(overall_proteins, picked_decoy_string_, picked_decoy_prefix_); }

    bool pepFDR = getStringOption_("FDR_type") == "PSM+peptide";
    // TODO Think about the implications of mixing PSMs from different files and searches.
    //   Score should be PEPs here. We could extract the original search scores, depending on preprocessing. PEPs allow some normalization but will
    //   disregard the absolute score differences between runs (i.e. if scores in one run are all lower than the ones in another run,
    //   do you want to filter them out preferably or do you say: this was a faulty run, if the decoys are equally bad, I want the
    //   best targets to be treated like the best targets from the other runs, even if the absolute match scores are much lower).
    if (pepFDR) { fdr.applyBasicPeptideLevel(consensus, true); }
    else { fdr.applyBasic(consensus, true); }

    if (! getFlag_("PeptideQuantification:quantify_decoys"))
    { // FDR filtering removed all decoy proteins -> update references and remove all unreferenced (decoy) PSMs
      IDFilter::removeDanglingProteinReferences(consensus, true);
      IDFilter::removeUnreferencedProteins(consensus, true); // if we don't filter peptides for now, we don't need this
      IDFilter::updateProteinGroups(overall_proteins.getIndistinguishableProteins(), overall_proteins.getHits());
      IDFilter::updateProteinGroups(overall_proteins.getProteinGroups(), overall_proteins.getHits());
    }

    // FDR filtering
    if (max_psm_fdr < 1.) // PSM level
    {
      for (auto& f : consensus)
      {
        IDFilter::filterHitsByScore(f.getPeptideIdentifications(), max_psm_fdr);
      }
      IDFilter::filterHitsByScore(consensus.getUnassignedPeptideIdentifications(), max_psm_fdr);
    }

    if (max_fdr < 1.) // protein level
    {
      IDFilter::filterHitsByScore(overall_proteins, max_fdr);
    }

    if (max_fdr < 1. || ! getFlag_("PeptideQuantification:quantify_decoys")) { IDFilter::removeDanglingProteinReferences(consensus, true); }

    if (max_psm_fdr < 1.) { IDFilter::removeUnreferencedProteins(consensus, true); }

    if (max_fdr < 1. || max_psm_fdr < 1. || ! getFlag_("PeptideQuantification:quantify_decoys"))
    {
      IDFilter::updateProteinGroups(overall_proteins.getIndistinguishableProteins(), overall_proteins.getHits());
      IDFilter::updateProteinGroups(overall_proteins.getProteinGroups(), overall_proteins.getHits());
    }

    if (overall_proteins.getHits().empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "No proteins left after FDR filtering. Please check the log and adjust your settings.");
    }

    // do we only want to keep strictly unique peptides (e.g., no groups)?
    // This filters for the VERY initially computed theoretical uniqueness calculated by PeptideIndexer
    //  which also means that e.g., target+decoy peptides are not unique
    if (! greedy_group_resolution && ! groups)
    {
      for (auto& f : consensus)
      {
        IDFilter::keepUniquePeptidesPerProtein(f.getPeptideIdentifications());
      }
      IDFilter::keepUniquePeptidesPerProtein(consensus.getUnassignedPeptideIdentifications());

      // Proteins whose peptides were all shared have no evidence left; drop them before grouping so
      // the groups below do not reference proteins that the cleanup after this function removes.
      IDFilter::removeUnreferencedProteins(consensus, true);

      if (overall_proteins.getHits().empty())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "No protein is supported by a strictly unique peptide. Note that "
                                            "'protein_quantification' = 'strictly_unique_peptides' requires the "
                                            "theoretical uniqueness annotated during peptide indexing (see '-fasta').");
      }

      // No indistinguishable groups were annotated above (that is what 'no groups' means here), but
      // quantification and every exporter report abundances on protein groups. Give each remaining
      // protein its own singleton group so a protein quantified from its strictly unique peptides is
      // actually reported instead of failing the "no indistinguishable protein groups" check later.
      overall_proteins.fillIndistinguishableGroupsWithSingletons();
    }

    // compute coverage (sequence was annotated during PeptideIndexing)
    // TODO: do you really want to compute coverage from unquantified peptides also?
    overall_proteins.computeCoverage(consensus, true);

    // TODO: this might not be correct if only the best peptidoform is kept
    // determine observed modifications (exclude fixed mods)
    overall_proteins.computeModifications(consensus, StringList(fixed_modifications.begin(), fixed_modifications.end()), true);

    return EXECUTION_OK;
  }


  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // Parameter handling
    //-------------------------------------------------------------

    // Read tool parameters
    const std::string out = getStringOption_("out");
    const std::string out_msstats = getStringOption_("out_msstats");
    const std::string out_cxml = getStringOption_("out_cxml");
    const std::string out_qpx = getOutputDirOption("out_qpx");
    feat_dir_ = getOutputDirOption("feat_dir");
    const bool detect_only = getFlag_("detect_only");

    // A detect-only invocation produces checkpoints and nothing else, so it has no result file to
    // require. Every other invocation still must be asked for something.
    if (! detect_only && out.empty() && out_msstats.empty() && out_cxml.empty() && out_qpx.empty())
    {
      throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                 "out/out_msstats/out_cxml/out_qpx");
    }
    if (detect_only && feat_dir_.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "'-detect_only' needs '-feat_dir' to write the checkpoints to.");
    }
    // Rejected rather than ignored: on its own the flag reads like "redo the work", and silently
    // doing nothing would look identical to a run that reused everything.
    if (getFlag_("force_recompute") && feat_dir_.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "'-force_recompute' only applies to '-feat_dir': without one "
                                        "there are no checkpoints to ignore, and every run is detected anyway.");
    }

    StringList in = getStringList_("in");
    StringList in_ids = getStringList_("ids");

    // '-in' and '-ids' are registered optional so that a combining invocation, which is covered
    // entirely by checkpoints, can omit them. Every other invocation still requires them, and
    // says so in the same terms it always did.
    if (feat_dir_.empty() && in.empty())
    {
      throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "in");
    }
    if (feat_dir_.empty() && in_ids.empty())
    {
      throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "ids");
    }

    std::string design_file = getStringOption_("design");
    std::string in_db = getStringOption_("fasta");

    if (! feat_dir_.empty())
    {
      // A run using checkpoints cannot fall back on a generated design: a detect-only invocation
      // sees one file and would generate Fraction_Group 1 for it, so every machine would
      // independently label its run the first one and the checkpoints would not compose.
      if (design_file.empty())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "'-feat_dir' requires an explicit '-design'. Without one the design is generated from '-in', "
          "which for a single-run invocation assigns Fraction_Group 1 to every run.");
      }
      // Without a FASTA there is no peptide indexing, hence no decoy affix, no theoretical
      // uniqueness and no protein sequences -- and a combining run has no way to obtain them later.
      if (in_db.empty())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "'-feat_dir' requires '-fasta': a checkpoint has to carry the peptide-indexing results, "
          "which a combining run cannot reconstruct.");
      }
      if (getStringOption_("quantification_method") == "spectral_counting")
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "'-feat_dir' does not apply to 'spectral_counting', which detects no features.");
      }
      if (! File::exists(feat_dir_) && ! File::makeDir(feat_dir_))
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, feat_dir_);
      }
      checkpoint_ctx_.fingerprint = checkpointFingerprint_();
      checkpoint_ctx_.fasta_stamp = inputStamp_(File::absolutePath(in_db));
      checkpoint_ctx_.openms_version = VersionInfo::getVersion();
      checkpoint_ctx_.openms_revision = VersionInfo::getRevision();
    }
    // Validate parameters
    if (in.size() != in_ids.size())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Number of spectra files (" + StringUtils::toStr(in.size()) + ") must match number of ID files (" + StringUtils::toStr(in_ids.size()) + ").");
    }
    if (getStringOption_("quantification_method") == "spectral_counting")
    {
      if (! out_msstats.empty())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "MSstats export for spectral counting data not supported. Please remove output file.");
      }
    }

    if (getStringOption_("targeted_only") != "false" && getStringOption_("pip_echo") != "false")
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "pip_echo requires targeted_only to be false");
    }

    if (getStringOption_("protein_quantification") == "strictly_unique_peptides" && in_db.empty())
    { // uniqueness comes from the 'protein_references' meta value, which is annotated by (re-)indexing
      OPENMS_LOG_WARN << "Warning: '-protein_quantification strictly_unique_peptides' filters peptides by their "
                         "theoretical uniqueness, which is annotated during peptide indexing. Without '-fasta' this "
                         "relies on the input identifications already being indexed (e.g., by PeptideIndexer); "
                         "otherwise no peptide will pass the filter." << endl;
    }

    //-------------------------------------------------------------
    // Experimental design: read or generate default
    //-------------------------------------------------------------
    if (! design_file.empty())
    { // load from file
      design_ = ExperimentalDesignFile::load(design_file, false);
    }
    else
    {
      OPENMS_LOG_INFO << "No experimental design file provided.\n"
                      << "Assuming a label-free experiment without fractionation.\n"
                      << endl;

      TextFile design_table;
      design_table.addLine("Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample\tMSstats_Condition\tMSstats_BioReplicate");

      Size count{1};
      for (std::string & s : in)
      {
        design_table.addLine(StringUtils::toStr(count) + "\t1\t" + s +"\t1\tSample" + StringUtils::toStr(count) + "\t" + StringUtils::toStr(count)+ "\t" + StringUtils::toStr(count));
        ++count;
      }
      design_ = ExperimentalDesignFile::load(design_table, false, "--no design file--");
    }

    // some sanity checks
    // extract basenames from experimental design and input files
    const auto& pl2fg = design_.getPathLabelToFractionGroupMapping(true);
    set<std::string> ed_basenames;
    for (const auto& p : pl2fg)
    {
      const std::string& filename = p.first.first;
      ed_basenames.insert(filename);
    }

    set<std::string> in_basenames;
    for (const auto & f : in)
    {
      const std::string& in_bn = File::basename(f);
      in_basenames.insert(in_bn);
    }

    if (! std::includes(ed_basenames.begin(), ed_basenames.end(), in_basenames.begin(), in_basenames.end()))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Spectra file basenames provided as input need to match a subset the experimental design file basenames.");
    }

    // With no '-in' the roster is the design itself: a combining run is covered entirely by
    // checkpoints, and filtering by an empty basename set would erase every row.
    Size nr_filtered = in_basenames.empty() ? 0 : design_.filterByBasenames(in_basenames);
    if (nr_filtered > 0)
    {
      OPENMS_LOG_WARN << "WARNING: " << nr_filtered
                      << " files from experimental design were not passed as mzMLs. Continuing with subset if the fractions still match."
                      << std::endl;
    }

    if (design_.getNumberOfLabels() != 1)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Experimental design is not label-free as it contains multiple labels.");
    }
    if (! design_.sameNrOfMSFilesPerFraction())
    {
      OPENMS_LOG_WARN << "WARNING: Different number of fractions for different samples provided. Support maybe limited in ProteomicsLFQ."
                      << std::endl;
    }

    std::map<unsigned int, std::vector<std::string> > frac2ms = design_.getFractionToMSFilesMapping();

    // experimental design file could contain URLs etc. that we want to overwrite with the actual input files
    for (auto& fraction_ms_files : frac2ms)
    {
      for (auto& s : fraction_ms_files.second)
      { // for all ms files of current fraction number
        // if basename in experimental design matches to basename in input file
        // overwrite experimental design to point to existing file (and only if they were different)
        if (auto it = std::find_if(in.begin(), in.end(),
              [&s] (const std::string& in_filename) { return File::basename(in_filename) == File::basename(s); }); // basename matches?
                 it != in.end() && s != *it) // and differ?
        {
          OPENMS_LOG_INFO << "Path of spectra files differ between experimental design (1) and input (2). Using the path of the input file as "
                          << "we know this file exists on the file system: '" << *it << "' vs. '" << s << endl;
          s = *it; // overwrite filename in design with filename in input files
        }
      }
    }

    for (auto& f : frac2ms)
    {
      writeDebug_("Fraction " + StringUtils::toStr(f.first) + ":", 10);
      for (const std::string & s : f.second)
      {
        writeDebug_("MS file: " + s, 10);
      }
    }

    // Map between mzML file and corresponding id file
    // We assume that these are provided in the exact same order.
    map<std::string, std::string> mzfile2idfile = DDAWorkflowCommons::mapMzML2Ids(in, in_ids);
    map<std::string, std::string> idfile2mzfile = DDAWorkflowCommons::mapId2MzMLs(mzfile2idfile);

    // TODO maybe check if mzMLs in experimental design match to mzMLs passed as in parameter
    //  IF both are present


    Param pep_param = getParam_().copy("Posterior Error Probability:", true);
    writeDebug_("Parameters passed to PEP algorithm", pep_param, 3);

    // TODO: inference parameter

    Param pq_param = getParam_().copy("ProteinQuantification:", true);
    writeDebug_("Parameters passed to PeptideAndProteinQuant algorithm", pq_param, 3);


    Param com_param = getParam_().copy("algorithm:common:", true);
    writeDebug_("Common parameters passed to both sub-algorithms (mtd and epd)", com_param, 3);

    set<std::string> fixed_modifications, variable_modifications;

    //-------------------------------------------------------------
    // Loading input
    //-------------------------------------------------------------
    ConsensusMap consensus;

    //-------------------------------------------------------------
    // feature-based quantifications
    //-------------------------------------------------------------
    if (getStringOption_("quantification_method") == "feature_intensity")
    {
      OPENMS_LOG_INFO << "Performing feature intensity-based quantification." << endl;
      for (auto const& ms_files : frac2ms) // for each fraction->ms file(s)
      {
        ConsensusMap consensus_fraction; // quantitative result for this fraction identifier

        ExitCodes e = quantifyFraction_(ms_files, mzfile2idfile, in_db, consensus_fraction, fixed_modifications, variable_modifications);

        if (e != EXECUTION_OK) { return e; }

        consensus.appendColumns(consensus_fraction); // append consensus map calculated for this fraction number
      } // end of scope of fraction related data

      if (detect_only)
      {
        OPENMS_LOG_INFO << "Feature checkpoints written to " << feat_dir_ << ": "
                        << runs_detected_ << " run(s) detected, " << checkpoints_reused_
                        << " reused. Stopping (-detect_only).\n";
        return EXECUTION_OK;
      }
      if (! feat_dir_.empty())
      {
        OPENMS_LOG_INFO << "Feature checkpoints: " << checkpoints_reused_ << " reused, "
                        << runs_detected_ << " detected.\n";
      }

      consensus.sortByPosition();
      consensus.sortPeptideIdentificationsByMapIndex();

      if (debug_level_ >= 666)
      {
        FileHandler().storeConsensusFeatures("debug_after_normalization.consensusXML", consensus, {FileTypes::CONSENSUSXML}, log_type_);
      }
    }
    else if (getStringOption_("quantification_method") == "spectral_counting")
    {
      OPENMS_LOG_INFO << "Performing spectral counting-based quantification." << endl;

      // init consensus map with basic experimental design information
      consensus.setExperimentType("label-free");

      auto& all_protein_ids = consensus.getProteinIdentifications();
      auto& all_peptide_ids = consensus.getUnassignedPeptideIdentifications();

      Size run_index(0);
      for (auto const& ms_files : frac2ms) // for each fraction->ms file(s) e.g.: Fraction1->FileA,FileB,FileC
      {
        const Size& fraction = ms_files.first;

        // debug output
        writeDebug_("Processing fraction number: " + StringUtils::toStr(fraction) + "\nFiles: ",  1);
        for (std::string const & mz_file : ms_files.second) { writeDebug_(mz_file,  1); }

        // for sanity checks we collect the primary MS run basenames as well as the ones stored in the ID files (below)
        StringList id_MS_run_ref;
        StringList in_MS_run = ms_files.second;
        const auto& path_label_to_fractiongroup = design_.getPathLabelToFractionGroupMapping(true);

        // for each MS file of current fraction (e.g., all MS files that measured the n-th fraction)
        for (std::string const & mz_file : ms_files.second)
        {
          const Size fraction_group = path_label_to_fractiongroup.at({File::basename(mz_file), 1});
          // load and clean identification data associated with MS run
          vector<ProteinIdentification> protein_ids;
          PeptideIdentificationList peptide_ids;
          const std::string& mz_file_abs_path = File::absolutePath(mz_file);
          const std::string& id_file_abs_path = File::absolutePath(mzfile2idfile.at(mz_file_abs_path));

          {
            std::string run_decoy_string;
            bool run_decoy_prefix = true;
            bool run_decoy_inferred = false;
            ExitCodes e = loadAndCleanupIDFile_(id_file_abs_path, mz_file, in_db, fraction_group, fraction, protein_ids, peptide_ids,
                                                fixed_modifications, variable_modifications,
                                                run_decoy_string, run_decoy_prefix, run_decoy_inferred);
            if (e != EXECUTION_OK) return e;
            recordDecoyAffix_(run_decoy_string, run_decoy_prefix, run_decoy_inferred, id_file_abs_path);
          }

          StringList id_msfile_ref;
          protein_ids[0].getPrimaryMSRunPath(id_msfile_ref);
          id_MS_run_ref.push_back(id_msfile_ref[0]);

          // append the ProteinIdentification run (contains backlink to MS file) and the PeptideIdentifications (PSMs for this fraction and MS run) to
          // the list of UnassignedPeptideIdentifications
          all_protein_ids.emplace_back(std::move(protein_ids[0]));
          all_peptide_ids.insert(all_peptide_ids.end(), std::make_move_iterator(peptide_ids.begin()), std::make_move_iterator(peptide_ids.end()));
        }

        ////////////////////////////////////////////////////////////
        // Annotate experimental design in consensus map
        ////////////////////////////////////////////////////////////
        // for each MS file (as provided in the experimental design)
        const auto& path_label_to_sampleidx = design_.getPathLabelToSampleMapping(true);
        for (std::string const & mz_file : ms_files.second)
        {
          const Size curr_fraction_group = path_label_to_fractiongroup.at({File::basename(mz_file), 1});
          consensus.getColumnHeaders()[run_index].label = "label-free";
          consensus.getColumnHeaders()[run_index].filename = mz_file;
          consensus.getColumnHeaders()[run_index].unique_id = 1 + run_index;
          consensus.getColumnHeaders()[run_index].setMetaValue("fraction", fraction);
          consensus.getColumnHeaders()[run_index].setMetaValue("fraction_group", curr_fraction_group);
          consensus.getColumnHeaders()[run_index].setMetaValue(
            "sample_name", design_.getSampleSection().getSampleName(path_label_to_sampleidx.at({File::basename(mz_file), 1})));
          ++run_index;
        }
      }
    }

    //-------------------------------------------------------------
    // ID related algorithms
    //-------------------------------------------------------------

    ExitCodes e = inferProteinGroups_(consensus, fixed_modifications);
    if (e != EXECUTION_OK) return e;

    // clean up references (assigned and unassigned)
    IDFilter::removeUnreferencedProteins(consensus, true);

    // only keep best scoring ID for each consensus feature
    IDConflictResolverAlgorithm::resolve(consensus);

    //-------------------------------------------------------------
    // Peptide quantification
    //-------------------------------------------------------------
    PeptideAndProteinQuant quantifier;

    // TODO Why is there no easy quantifier.run(consensus,[inference_prot_ids]) function??
    if (getStringOption_("quantification_method") == "feature_intensity")
    {
      quantifier.setParameters(pq_param);
      quantifier.readQuantData(consensus, design_);
    }
    else if (getStringOption_("quantification_method") == "spectral_counting")
    {
      pq_param.setValue("top:aggregate", "sum");
      pq_param.setValue("top:N", 0); // all
      pq_param.setValue("consensus:normalize", "false");
      quantifier.setParameters(pq_param);

      quantifier.readQuantData(consensus.getProteinIdentifications(), consensus.getUnassignedPeptideIdentifications(), design_);
    }

    // nothing to filter. everything in consensus should be uptodate with inference.
    // on peptide level it does not annotate anything anyway
    quantifier.quantifyPeptides();

    //-------------------------------------------------------------
    // Protein quantification
    //-------------------------------------------------------------

    // Should always be there by now, even if just singletons
    ProteinIdentification& inferred_proteins = consensus.getProteinIdentifications()[0];
    if (inferred_proteins.getIndistinguishableProteins().empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No information on indistinguishable protein groups found.");
    }

    quantifier.quantifyProteins(inferred_proteins);
    auto const& protein_quants = quantifier.getProteinResults();
    if (protein_quants.empty()) { OPENMS_LOG_WARN << "Warning: No proteins were quantified." << endl; }

    //-------------------------------------------------------------
    // Export of MzTab file as final output
    //-------------------------------------------------------------

    // Annotate quants to protein(groups) for easier export in mzTab
    // Note: we keep protein groups that have not been quantified
    quantifier.annotateQuantificationsToProteins(protein_quants, inferred_proteins, false);

    // For correctness, we would need to set the run reference in the pepIDs of the consensusXML all to the first run then
    // And probably make sure that peptides that correspond to filtered out proteins are not producing errors
    // e.g. by removing them with a Filter beforehand.

    consensus.resolveUniqueIdConflicts(); // TODO: find out if this is still needed

    {
      if (!out_qpx.empty())
      {
        OPENMS_LOG_INFO << "Exporting QPX Parquet files to: " << out_qpx << std::endl;

        // Validate the whole collection before the first write. Each view refuses input it
        // cannot represent, but only before its OWN file, so a refusal in the psm or pg view
        // would leave the earlier files behind as a partial collection.
        if (!QPXCollectionExport::requireExportable(consensus, design_))
        {
          return CANNOT_WRITE_OUTPUT_FILE; // already logged
        }

        // Whatever the preflight cannot decide up front - a row-level refusal, an I/O error -
        // is undone here: without the commit below, every file already written is removed.
        QPXCollectionExport::Transaction qpx_collection(out_qpx);

        // Collected while the feature rows are built, then handed to the psm view so it can fill
        // psm.feature_id. One pass produces both directions, which is what keeps them reciprocal.
        QPXIdentity::FeatureLinks feature_links;

        // Feature-level export
        if (!ConsensusMapArrowExport::exportToParquet(consensus, out_qpx + "/quantms.feature.parquet",
                                                      ParquetWriteConfig{}, &feature_links))
        {
          OPENMS_LOG_ERROR << "Failed to write features Parquet file" << std::endl;
          return CANNOT_WRITE_OUTPUT_FILE;
        }

        // PSM-level export: collect all peptide IDs from consensus map
        PeptideIdentificationList all_pepids;
        for (const auto& feature : consensus)
        {
          for (const auto& pepid : feature.getPeptideIdentifications())
          {
            all_pepids.push_back(pepid);
          }
        }
        for (const auto& pepid : consensus.getUnassignedPeptideIdentifications())
        {
          all_pepids.push_back(pepid);
        }

        if (!QPXFile::exportToParquet(consensus.getProteinIdentifications(), all_pepids,
                                      out_qpx + "/quantms.psm.parquet", /*export_all_psms=*/false,
                                      ParquetWriteConfig{}, &feature_links))
        {
          OPENMS_LOG_ERROR << "Failed to write PSM Parquet file" << std::endl;
          return CANNOT_WRITE_OUTPUT_FILE;
        }

        // Protein group export
        // Pass the design that drove quantification: QPX 1.1 keys the pg view on the set of
        // files aggregated into one quantity, and the design is what defines that grouping and
        // the sample numbering the protein abundances are stored under.
        if (!ProteinGroupArrowExport::exportToParquet(consensus, design_, out_qpx + "/quantms.pg.parquet"))
        {
          OPENMS_LOG_ERROR << "Failed to write protein groups Parquet file" << std::endl;
          return CANNOT_WRITE_OUTPUT_FILE;
        }

        qpx_collection.commit(); // all three views written
      }
    }

    if (!out_cxml.empty())
    {
      // Note: idXML and consensusXML doesn't support writing quantification at protein groups
      // (they are nevertheless stored and passed to mzTab for proper export)
      FileHandler().storeConsensusFeatures(out_cxml, consensus, {FileTypes::CONSENSUSXML}, log_type_);
    }

    if (!out.empty())
    {
      // Fill mzTab with meta data and quants annotated in identification data structure
      const bool report_unidentified_features(false);
      const bool report_unmapped(true); // TODO we should make a distinction from unassigned after conflict resolution and unassigned because unmappable
      const bool report_subfeatures(false);
      const bool report_unidentified_spectra(false);
      const bool report_not_only_best_psm_per_spectrum(false);

      MzTabFile().store(out, consensus,
                        false, // first run is inference but also a properly merged run, so we don't need the hack
                        report_unidentified_features, report_unmapped, report_subfeatures, report_unidentified_spectra,
                        report_not_only_best_psm_per_spectrum);
    }

    if (! out_msstats.empty())
    {
      IDFilter::removeEmptyIdentifications(consensus); // MzTab stream exporter currently doesn't support IDs with empty hits.

      MSstatsFile msstats;
      // TODO: add a helper method to quickly check if experimental design file contain the right columns
      //  (and put this at start of tool)

      // shrink protein runs to the one containing the inference data
      consensus.getProteinIdentifications().resize(1);

      msstats.storeLFQ(out_msstats, consensus, design_, StringList(),
                       false, // lfq
                       "MSstats_BioReplicate", "MSstats_Condition", "max");
    }


    return EXECUTION_OK;
  }
  std::string picked_decoy_string_;
  bool picked_decoy_prefix_ = true;
  bool picked_decoy_seen_ = false; ///< has an affix been inferred from an input file yet?
  std::string feat_dir_;                  ///< empty unless -feat_dir was given
  CheckpointContext checkpoint_ctx_;      ///< the run-independent half, built once in main_
  Size checkpoints_reused_ = 0;
  Size runs_detected_ = 0;
  ExperimentalDesign design_;
};

int main(int argc, const char** argv)
{
  ProteomicsLFQ tool;
  return tool.main(argc, argv);
}

/// @endcond
