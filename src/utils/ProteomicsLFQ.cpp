// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/FORMAT/MSstatsFile.h>

#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>
#include <OpenMS/DATASTRUCTURES/CalibrationData.h>
#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/FILTERING/CALIBRATION/MZTrafoModel.h>
#include <OpenMS/FILTERING/CALIBRATION/PrecursorCorrection.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h>
#include <OpenMS/FILTERING/DATAREDUCTION/FeatureFindingMetabo.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
//#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
#include <OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/BayesianProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDBoostGraph.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/PeptideProteinResolution.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h>
#include <OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h>


#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>

#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/METADATA/SpectrumMetaDataLookup.h>

#include <OpenMS/KERNEL/ConversionHelper.h>

#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>

#include <OpenMS/FILTERING/ID/IDFilter.h>

#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h>

using namespace OpenMS;
using namespace std;
using Internal::IDBoostGraph;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_ProteomicsLFQ ProteomicsLFQ

  ProteomicsLFQ performs label-free quantification of peptides and proteins. @n

  Input: @n
         - Spectra in mzML format
         - Identifications in idXML or mzIdentML format with posterior error probabilities
           as score type.
           To generate those we suggest to run:
           1. PeptideIndexer to annotate target and decoy information.
           2. PSMFeatureExtractor to annotate percolator features.
           3. PercolatorAdapter tool (score_type = 'q-value', -post-processing-tdc)
           4. IDFilter (pep:score = 0.01) to filter PSMs at 1% FDR
           
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
            Transfer of Ids (match between runs) is performed by transfering feature identifications to coeluting, unidentified features with similar mass and RT in other runs.

         Requantification:
         2. Optionally, a requantification step is performed that tries to fill NA values.
            If a peptide has been quantified in more than half of all maps, the peptide is selected for requantification.
            In that case, the mean observed RT (and theoretical m/z) of the peptide is used to perform a second round of targeted extraction.
  Output:
         - mzTab file with analysis results
         - MSstats file with analysis results for statistical downstream analysis in MSstats
         - ConsensusXML file for visualization and further processing in OpenMS

  // experiments TODO:
  // - change percentage of missingness in ID transfer
  // - disable elution peak fit

  Potential scripts to perform the search can be found under src/tests/topp/ProteomicsLFQTestScripts
 **/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class ProteomicsLFQ :
  public TOPPBase
{
public:
  ProteomicsLFQ() :
    TOPPBase("ProteomicsLFQ", "A standard proteomics LFQ pipeline.", false)
  {
  }

protected:  
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<file list>", StringList(), "Input files");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerInputFileList_("ids", "<file list>", StringList(), 
      "Identifications filtered at PSM level (e.g., q-value < 0.01)."
      "And annotated with PEP as main score.\n"
      "We suggest using:\n"
      "1. PeptideIndexer to annotate target and decoy information.\n"
      "2. PSMFeatureExtractor to annotate percolator features.\n"
      "3. PercolatorAdapter tool (score_type = 'q-value', -post-processing-tdc)\n"
      "4. IDFilter (pep:score = 0.01)\n"
      "To obtain well calibrated PEPs and an inital reduction of PSMs\n"
      "ID files must be provided in same order as spectra files.");
    setValidFormats_("ids", ListUtils::create<String>("idXML,mzId"));

    registerInputFile_("design", "<file>", "", "design file", false);
    setValidFormats_("design", ListUtils::create<String>("tsv"));

    registerInputFile_("fasta", "<file>", "", "fasta file", false);
    setValidFormats_("fasta", ListUtils::create<String>("fasta"));

    registerOutputFile_("out", "<file>", "", "output mzTab file");
    setValidFormats_("out", ListUtils::create<String>("mzTab"));

    registerOutputFile_("out_msstats", "<file>", "", "output MSstats input file", false, false);
    setValidFormats_("out_msstats", ListUtils::create<String>("csv"));

    registerOutputFile_("out_cxml", "<file>", "", "output consensusXML file", false, false);
    setValidFormats_("out_cxml", ListUtils::create<String>("consensusXML"));

    registerDoubleOption_("proteinFDR", "<threshold>", 0.05, "Protein FDR threshold (0.05=5%).", false);
    setMinFloat_("proteinFDR", 0.0);
    setMaxFloat_("proteinFDR", 1.0);

    registerDoubleOption_("seedThreshold", "<threshold>", 1e4, "Peak intensity threshold applied in seed detection.", false, true);

    registerDoubleOption_("psmFDR", "<threshold>", 1.0, "PSM FDR threshold (e.g. 0.05=5%)."
                          " If Bayesian inference was chosen, it is equivalent with a peptide FDR", false);
    setMinFloat_("psmFDR", 0.0);
    setMaxFloat_("psmFDR", 1.0);

    //TODO expose all parameters of the inference algorithms (e.g. aggregation methods etc.)?
    registerStringOption_("protein_inference", "<option>", "aggregation",
      "Infer proteins:\n" 
      "aggregation  = aggregates all peptide scores across a protein (by calculating the maximum) \n"
      "bayesian     = computes a posterior probability for every protein based on a Bayesian network.\n"
      "               Note: 'bayesian' only uses and reports the best PSM per peptide.",
      false, true);
    setValidStrings_("protein_inference", ListUtils::create<String>("aggregation,bayesian"));

    registerStringOption_("protein_quantification", "<option>", "unique_peptides",
      "Quantify proteins based on:\n"
      "unique_peptides = use peptides mapping to single proteins or a group of indistinguishable proteins"
      "(according to the set of experimentally identified peptides).\n"
      "strictly_unique_peptides = use peptides mapping to a unique single protein only.\n"
      "shared_peptides = use shared peptides only for its best group (by inference score)", false, true);
    setValidStrings_("protein_quantification", ListUtils::create<String>("unique_peptides,strictly_unique_peptides,shared_peptides"));
    registerStringOption_("quantification_method", "<option>", 
      "feature_intensity", 
      "feature_intensity: MS1 signal.\n"
      "spectral_counting: PSM counts.", false, false);
    setValidStrings_("quantification_method", ListUtils::create<String>("feature_intensity,spectral_counting"));

    registerStringOption_("targeted_only", "<option>", "false", 
      "true: Only ID based quantification.\n"
      "false: include unidentified features so they can be linked to identified ones (=match between runs).", false, false);
    setValidStrings_("targeted_only", ListUtils::create<String>("true,false"));

    // TODO: support transfer with SVM if we figure out a computational efficient way to do it.
    registerStringOption_("transfer_ids", "<option>", "false", 
      "Requantification using mean of aligned RTs of a peptide feature.\n"
      "Only applies to peptides that were quantified in more than 50% of all runs (of a fraction).", false, false);
    setValidStrings_("transfer_ids", ListUtils::create<String>("false,mean"));

    registerStringOption_("mass_recalibration", "<option>", "false", "Mass recalibration.", false, true);
    setValidStrings_("mass_recalibration", ListUtils::create<String>("true,false"));


    /// TODO: think about export of quality control files (qcML?)

    Param pp_defaults = PeakPickerHiRes().getDefaults();
    for (const auto& s : {"report_FWHM", "report_FWHM_unit", "SignalToNoise:win_len", "SignalToNoise:bin_count", "SignalToNoise:min_required_elements", "SignalToNoise:write_log_messages"} )
    {
      pp_defaults.addTag(s, "advanced");
    }

    Param ffi_defaults = FeatureFinderIdentificationAlgorithm().getDefaults();
    ffi_defaults.setValue("svm:samples", 10000); // restrict number of samples for training
    ffi_defaults.setValue("svm:log2_C", DoubleList({-2.0, 5.0, 15.0})); 
    ffi_defaults.setValue("svm:log2_gamma", DoubleList({-3.0, -1.0, 2.0})); 
    ffi_defaults.setValue("svm:min_prob", 0.9); // keep only feature candidates with > 0.9 probability of correctness
    // hide entries
    for (const auto& s : {"svm:samples", "svm:log2_C", "svm:log2_gamma", "svm:min_prob", "svm:no_selection", "svm:xval_out", "svm:kernel", "svm:xval", "candidates_out", "extract:n_isotopes", "model:type"} )
    {
      ffi_defaults.addTag(s, "advanced");
    }
    ffi_defaults.remove("detect:peak_width"); // set from data

    Param ma_defaults = MapAlignmentAlgorithmIdentification().getDefaults();
    ma_defaults.setValue("max_rt_shift", 0.1);
    ma_defaults.setValue("use_unassigned_peptides", "false");
    ma_defaults.setValue("use_feature_rt", "true");

    // hide entries
    for (const auto& s : {"use_unassigned_peptides", "use_feature_rt", "score_cutoff", "min_score"} )
    {
      ma_defaults.addTag(s, "advanced");
    }

    //Param fl_defaults = FeatureGroupingAlgorithmKD().getDefaults();
    Param fl_defaults = FeatureGroupingAlgorithmQT().getDefaults();
    fl_defaults.setValue("distance_MZ:max_difference", 10.0);
    fl_defaults.setValue("distance_MZ:unit", "ppm");
    fl_defaults.setValue("distance_MZ:weight", 5.0);
    fl_defaults.setValue("distance_intensity:weight", 0.1); 
    fl_defaults.setValue("use_identifications", "true"); 
    fl_defaults.remove("distance_RT:max_difference"); // estimated from data
    for (const auto& s : {"distance_MZ:weight", "distance_intensity:weight", "use_identifications", "ignore_charge", "ignore_adduct"} )
    {
      fl_defaults.addTag(s, "advanced");
    }

    Param pq_defaults = PeptideAndProteinQuant().getDefaults();
    // overwrite algorithm default so we export everything (important for copying back MSstats results)
    pq_defaults.setValue("include_all", "true"); 
    pq_defaults.addTag("include_all", "advanced");

    // combine parameters of the individual algorithms
    Param combined;
    combined.insert("Centroiding:", pp_defaults);
    combined.insert("PeptideQuantification:", ffi_defaults);
    combined.insert("Alignment:", ma_defaults);
    combined.insert("Linking:", fl_defaults);
    combined.insert("ProteinQuantification:", pq_defaults);

    registerFullParam_(combined);
  }

  // Map between mzML file and corresponding id file
  // Warn if the primaryMSRun indicates that files were provided in the wrong order.
  map<String, String> mapMzML2Ids_(StringList & in, StringList & in_ids)
  {
    // Detect the common case that ID files have same names as spectra files
    if (!File::validateMatchingFileNames(in, in_ids, true, true, false)) // only basenames, without extension, only order
    {
      // Spectra and id files have the same set of basenames but appear in different order. -> this is most likely an error
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "ID and spectra file match but order of file names seem to differ. They need to be provided in the same order.");
    }

    map<String, String> mzfile2idfile;
    for (Size i = 0; i != in.size(); ++i)
    {
      const String& in_abs_path = File::absolutePath(in[i]);
      const String& id_abs_path = File::absolutePath(in_ids[i]);
      mzfile2idfile[in_abs_path] = id_abs_path;      
      writeDebug_("Spectra: " + in[i] + "\t Ids: " + in_ids[i],  1);
    }
    return mzfile2idfile;
  }

  // map back
  map<String, String> mapId2MzMLs_(const map<String, String>& m2i)
  {
    map<String, String> idfile2mzfile;
    for (const auto& m : m2i)
    {
      idfile2mzfile[m.second] = m.first;
    }
    return idfile2mzfile;
  }

  ExitCodes centroidAndCorrectPrecursors_(const String & mz_file, MSExperiment & ms_centroided)
  { 
    Param pp_param = getParam_().copy("Centroiding:", true);
    writeDebug_("Parameters passed to PeakPickerHiRes algorithm", pp_param, 3);

    // create scope for raw data so it is properly freed (Note: clear() is not sufficient)
    // load raw file
    MzMLFile mzML_file;
    mzML_file.setLogType(log_type_);

    PeakMap ms_raw;
    mzML_file.load(mz_file, ms_raw);
    ms_raw.clearMetaDataArrays();

    if (ms_raw.empty())
    {
      OPENMS_LOG_WARN << "The given file does not contain any spectra.";
      return INCOMPATIBLE_INPUT_DATA;
    }

    // remove MS2 peak data and check if spectra are sorted
    for (Size i = 0; i < ms_raw.size(); ++i)
    {
      if (ms_raw[i].getMSLevel() == 2)
      {
        ms_raw[i].clear(false);  // delete MS2 peaks
      }
      if (!ms_raw[i].isSorted())
      {
        ms_raw[i].sortByPosition();
        writeLog_("Info: Sorted peaks by m/z.");
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
    writeLog_("Info: Corrected " + String(corrected_to_highest_intensity_peak.size()) + " precursors.");
    if (!deltaMZs.empty())
    {
      vector<double> deltaMZs_ppm, deltaMZs_ppmabs;
      for (Size i = 0; i != deltaMZs.size(); ++i)
      {
        deltaMZs_ppm.push_back(Math::getPPM(mzs[i], mzs[i] + deltaMZs[i]));
        deltaMZs_ppmabs.push_back(Math::getPPMAbs(mzs[i], mzs[i] + deltaMZs[i]));
      }

      double median = Math::median(deltaMZs_ppm.begin(), deltaMZs_ppm.end());
      double MAD =  Math::MAD(deltaMZs_ppm.begin(), deltaMZs_ppm.end(), median);
      double median_abs = Math::median(deltaMZs_ppmabs.begin(), deltaMZs_ppmabs.end());
      double MAD_abs = Math::MAD(deltaMZs_ppmabs.begin(), deltaMZs_ppmabs.end(), median_abs);
      writeLog_("Precursor correction:\n  median        = " 
        + String(median) + " ppm  MAD = " + String(MAD)
        + "\n  median (abs.) = " + String(median_abs) 
        + " ppm  MAD = " + String(MAD_abs));
    }
    return EXECUTION_OK;
  }

  void recalibrateMasses_(MSExperiment & ms_centroided, vector<PeptideIdentification>& peptide_ids, const String & id_file_abs_path)
  {
    InternalCalibration ic;
    ic.setLogType(log_type_);
    ic.fillCalibrants(peptide_ids, 25.0); // >25 ppm maximum deviation defines an outlier TODO: check if we need to adapt this
    if (ic.getCalibrationPoints().size() <= 1) return;

    // choose calibration model based on number of calibration points

    // there seem to be some problems with the QUADRATIC model that we first need to investigate
    //MZTrafoModel::MODELTYPE md = (ic.getCalibrationPoints().size() == 2) ? MZTrafoModel::LINEAR : MZTrafoModel::QUADRATIC;
    //bool use_RANSAC = (md == MZTrafoModel::LINEAR || md == MZTrafoModel::QUADRATIC);

    MZTrafoModel::MODELTYPE md = MZTrafoModel::LINEAR;
    bool use_RANSAC = true;

    Size RANSAC_initial_points = (md == MZTrafoModel::LINEAR) ? 2 : 3;
    Math::RANSACParam p(RANSAC_initial_points, 70, 10, 30, true); // TODO: check defaults (taken from tool)
    MZTrafoModel::setRANSACParams(p);
    // these limits are a little loose, but should prevent grossly wrong models without burdening the user with yet another parameter.
    MZTrafoModel::setCoefficientLimits(25.0, 25.0, 0.5); 

    IntList ms_level = {1};
    double rt_chunk = 300.0; // 5 minutes
    String qc_residual_path, qc_residual_png_path;
    if (debug_level_ >= 1)
    {
      const String & id_basename = File::basename(id_file_abs_path);
      qc_residual_path = id_basename + "qc_residuals.tsv";
      qc_residual_png_path = id_basename + "qc_residuals.png";
    } 

    if (!ic.calibrate(ms_centroided, ms_level, md, rt_chunk, use_RANSAC, 
                  10.0,
                  5.0, 
                  "",                      
                  "",
                  qc_residual_path,
                  qc_residual_png_path,
                  "Rscript"))
    {
      OPENMS_LOG_WARN << "\nCalibration failed. See error message above!" << std::endl;
    }
  }

  double estimateMedianChromatographicFWHM_(MSExperiment & ms_centroided)
  {
    MassTraceDetection mt_ext;
    Param mtd_param = mt_ext.getParameters();
    writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);

    std::vector<MassTrace> m_traces;
    mt_ext.run(ms_centroided, m_traces, 1000);

    std::vector<double> fwhm_1000;
    for (auto &m : m_traces)
    {
      if (m.getSize() == 0) continue;
      m.updateMeanMZ();
      m.updateWeightedMZsd();
      double fwhm = m.estimateFWHM(false);
      fwhm_1000.push_back(fwhm);
    }

    double median_fwhm = Math::median(fwhm_1000.begin(), fwhm_1000.end());

    OPENMS_LOG_INFO << "Median chromatographic FWHM: " << median_fwhm << std::endl;

    return median_fwhm;
  }

  void calculateSeeds_(const MSExperiment & ms_centroided, FeatureMap & seeds, double median_fwhm)
  {
    MSExperiment e;
    for (const auto& s : ms_centroided)
    { 
      if (s.getMSLevel() == 1) 
      {              
        e.addSpectrum(s);
      }
    }
    ThresholdMower threshold_mower_filter;
    Param tm = threshold_mower_filter.getParameters();
    tm.setValue("threshold", getDoubleOption_("seedThreshold"));  // TODO: derive from data
    threshold_mower_filter.setParameters(tm);
    threshold_mower_filter.filterPeakMap(e);

    FeatureFinderMultiplexAlgorithm algorithm;
    Param p = algorithm.getParameters();
    p.setValue("algorithm:labels", "");
    p.setValue("algorithm:charge", "2:5");
    p.setValue("algorithm:rt_typical", median_fwhm * 3.0);
    p.setValue("algorithm:rt_band", 3.0); // max 3 seconds shifts between isotopic traces
    p.setValue("algorithm:rt_min", median_fwhm * 0.5);
    p.setValue("algorithm:spectrum_type", "centroid");
    algorithm.setParameters(p);
    const bool progress(true);
    algorithm.run(e, progress);
    seeds = algorithm.getFeatureMap(); 
    OPENMS_LOG_INFO << "Using " << seeds.size() << " seeds from untargeted feature extraction." << endl;
  }


  // aligns the feature maps 
  double align_(
    vector<FeatureMap> & feature_maps, 
    vector<TransformationDescription>& transformations
  )
  {
    if (feature_maps.size() > 1) // do we have several maps to align / link?
    {
      Param ma_param = getParam_().copy("Alignment:", true);
      writeDebug_("Parameters passed to MapAlignmentAlgorithmIdentification algorithm", ma_param, 3);

      // Determine reference from data, otherwise a change in order of input files
      // leads to slightly different results
      const int reference_index(-1); // set no reference (determine from data)
      MapAlignmentAlgorithmIdentification aligner;
      aligner.setLogType(log_type_);
      aligner.setParameters(ma_param);

      Param model_params = TOPPMapAlignerBase::getModelDefaults("b_spline");
      String model_type = model_params.getValue("type");
      model_params = model_params.copy(model_type + ":", true);

      try
      {
        aligner.align(feature_maps, transformations, reference_index);
      }
      catch (Exception::MissingInformation& err)
      {
        if (getFlag_("force"))
        {
          OPENMS_LOG_ERROR
            << "Error: alignment failed. Details:\n" << err.getMessage()
            << "\nProcessing will continue using 'identity' transformations."
            << endl;
          model_type = "identity";
          transformations.resize(feature_maps.size());
        }
        else throw;
      }

      // find model parameters (if model_type == "identity" the fit is a NOP):
      vector<TransformationDescription::TransformationStatistics> alignment_stats;
      for (TransformationDescription & t : transformations)
      {
        writeDebug_("Using " + String(t.getDataPoints().size()) + " points in fit.", 1); 
        if (t.getDataPoints().size() > 10)
        {
          t.fitModel(model_type, model_params);
        }
        t.printSummary(OpenMS_Log_debug);
        alignment_stats.emplace_back(t.getStatistics());
      }

      // determine maximum RT shift after transformation that includes all high confidence IDs 
      using TrafoStat = TransformationDescription::TransformationStatistics;
      for (auto & s : alignment_stats)
      {
        OPENMS_LOG_INFO << "Alignment differences (second) for percentiles (before & after): " << endl;
        OPENMS_LOG_INFO << ListUtils::concatenate(s.percents,"%\t") << "%" << endl;
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

      double max_alignment_diff = std::max_element(alignment_stats.begin(), alignment_stats.end(),
              [](TrafoStat a, TrafoStat b) 
              { return a.percentiles_after[100] > b.percentiles_after[100]; })->percentiles_after[100];
      // sometimes, very good alignments might lead to bad overall performance. Choose 2 minutes as minimum.
      OPENMS_LOG_INFO << "Max alignment difference (seconds): " << max_alignment_diff << endl;
      max_alignment_diff = std::max(max_alignment_diff, 120.0);
      return max_alignment_diff;
    }
    return 0;
  }

  void transform_(
    vector<FeatureMap>& feature_maps, 
    vector<TransformationDescription>& transformations
  )
  {
    if (feature_maps.size() > 1 && !transformations.empty())
    {
      // Apply transformations
      for (Size i = 0; i < feature_maps.size(); ++i)
      {
        try 
        {
          MapAlignmentTransformer::transformRetentionTimes(feature_maps[i],
            transformations[i]);
        } catch (Exception::IllegalArgument& e)
        {
          OPENMS_LOG_WARN << e.getMessage() << endl;
        }
          
        if (debug_level_ > 666)
        {
          // plot with e.g.:
          // Rscript ../share/OpenMS/SCRIPTS/plot_trafo.R debug_trafo_1.trafoXML debug_trafo_1.pdf
          TransformationXMLFile().store("debug_trafo_" + String(i) + ".trafoXML", transformations[i]);
        }
      }
    }
  }

  //-------------------------------------------------------------
  // Link all features of this fraction
  //-------------------------------------------------------------
  void link_(
    vector<FeatureMap> & feature_maps, 
    double median_fwhm,
    double max_alignment_diff,
    ConsensusMap & consensus_fraction
  )
  {
    Param fl_param = getParam_().copy("Linking:", true);
    writeDebug_("Parameters passed to feature grouping algorithm", fl_param, 3);

    writeDebug_("Linking: " + String(feature_maps.size()) + " features.", 1);

    // grouping tolerance = max alignment error + median FWHM
    FeatureGroupingAlgorithmQT linker;
    fl_param.setValue("distance_RT:max_difference", 2.0 * max_alignment_diff + 2.0 * median_fwhm);
    linker.setParameters(fl_param);      
/*
    FeatureGroupingAlgorithmKD linker;
    fl_param.setValue("warp:rt_tol", 2.0 * max_alignment_diff + 2.0 * median_fwhm);
    fl_param.setValue("link:rt_tol", 2.0 * max_alignment_diff + 2.0 * median_fwhm);
    fl_param.setValue("link:mz_tol", 10.0);
    fl_param.setValue("mz_unit", "ppm");
    linker.setParameters(fl_param);      
*/
    linker.group(feature_maps, consensus_fraction);
    OPENMS_LOG_INFO << "Size of consensus fraction: " << consensus_fraction.size() << endl;
    assert(!consensus_fraction.empty());
  }

  // Align and link.
  // @return maximum alignment difference observed (to guide linking)
  double alignAndLink_(
    vector<FeatureMap> & feature_maps, 
    ConsensusMap & consensus_fraction,
    vector<TransformationDescription>& transformations,
    const double median_fwhm)
  {
    double max_alignment_diff(0.0);

    if (feature_maps.size() > 1)
    {
      max_alignment_diff = align_(feature_maps, transformations);

      transform_(feature_maps, transformations);

      link_(feature_maps, 
        median_fwhm, 
        max_alignment_diff, 
        consensus_fraction);
    }
    else // only one feature map
    {
      MapConversion::convert(0, feature_maps.back(), consensus_fraction);                           
    }

    return max_alignment_diff;
  }

  // determine cooccurance of peptide in different runs
  // returns map sequence+charge -> map index in consensus map 
  map<pair<String, UInt>, vector<int> > getPeptideOccurrence_(const ConsensusMap &cons)
  {
    map<Size, UInt> num_consfeat_of_size;
    map<Size, UInt> num_consfeat_of_size_with_id;

    map<pair<String, UInt>, vector<int> > seq_charge2map_occurence;
    for (ConsensusMap::const_iterator cmit = cons.begin(); cmit != cons.end(); ++cmit)
    {
      ++num_consfeat_of_size[cmit->size()];
      const auto& pids = cmit->getPeptideIdentifications();
      if (!pids.empty())
      {
        ++num_consfeat_of_size_with_id[cmit->size()];

        // count how often a peptide/charge pair has been observed in the different maps
        const vector<PeptideHit>& phits = pids[0].getHits();
        if (!phits.empty())
        {
          const String s = phits[0].getSequence().toString();
          const int z = phits[0].getCharge();

          if (seq_charge2map_occurence[make_pair(s,z)].empty())
          {
            seq_charge2map_occurence[make_pair(s,z)] = vector<int>(cons.getColumnHeaders().size(), 0);
          }

          // assign id to all dimensions in the consensus feature
          for (auto const & f : cmit->getFeatures())
          {
            Size map_index = f.getMapIndex();
            seq_charge2map_occurence[make_pair(s,z)][map_index] += 1;
          }
        }
      }
    }
    return seq_charge2map_occurence;
  }

  // simple transfer between runs
  // if a peptide has not been quantified in more than min_occurrence runs, then take all consensus features that have it identified at least once
  // and transfer the ID with RT of the the consensus feature (the average if we have multiple consensus elements)
  multimap<Size, PeptideIdentification> transferIDsBetweenSameFraction_(const ConsensusMap& consensus_fraction, Size min_occurrence = 3)
  {
    // determine occurrence of ids
    map<pair<String, UInt>, vector<int> > occurrence = getPeptideOccurrence_(consensus_fraction);

    // build map of missing ids
    map<pair<String, UInt>, set<int> > missing; // set of maps missing the id
    for (auto & o : occurrence)
    {
      // more than min_occurrence elements in consensus map that are non-zero?
      const Size count_non_zero = (Size) std::count_if(o.second.begin(), o.second.end(), [](int i){return i > 0;});

      if (count_non_zero >= min_occurrence
       && count_non_zero < o.second.size())
      {
        for (Size i = 0; i != o.second.size(); ++i)
        {
          // missing ID for this consensus element
          if (o.second[i] == 0) { missing[o.first].insert(i); }
        }
      }
    }

    Size n_transferred_ids(0);
    // create representative id to transfer to missing
    multimap<Size, PeptideIdentification> transfer_ids;
    for (auto & c : consensus_fraction)
    {
      const auto& pids = c.getPeptideIdentifications();
      if (pids.empty()) continue; // skip consensus feature without IDs 

      const vector<PeptideHit>& phits = pids[0].getHits();
      if (phits.empty()) continue; // skip no PSM annotated

      const String s = phits[0].getSequence().toString();
      const int z = phits[0].getCharge();
      pair<String, UInt> seq_z = make_pair(s, z);
      map<pair<String, UInt>, set<int> >::const_iterator it = missing.find(seq_z);

      if (it == missing.end()) continue; // skip sequence and charge not marked as missing in one of the other maps

      for (int idx : it->second)
      {
        // use consensus feature ID and retention time to transfer between runs
        pair<Size, PeptideIdentification> p = make_pair(idx, pids[0]);
        p.second.setRT(c.getRT());
        transfer_ids.insert(p);
        ++n_transferred_ids;
      }
    }
    OPENMS_LOG_INFO << "Transfered IDs: " << n_transferred_ids << endl;
    return transfer_ids;
  }

  ExitCodes checkSingleRunPerID_(const vector<ProteinIdentification>& protein_ids, const String& id_file_abs_path)
  {
    if (protein_ids.size() != 1)
    {
      OPENMS_LOG_FATAL_ERROR << "Exactly one protein identification run must be annotated in " << id_file_abs_path << endl;
      return ExitCodes::INCOMPATIBLE_INPUT_DATA;
    }

    StringList run_paths;
    protein_ids[0].getPrimaryMSRunPath(run_paths);
    if (run_paths.size() != 1)
    {
      OPENMS_LOG_FATAL_ERROR << "ProteomicsLFQ does not support merged ID runs. ID file: " << id_file_abs_path << endl;
      return ExitCodes::INCOMPATIBLE_INPUT_DATA;
    }
    
    return EXECUTION_OK;
  }

  ExitCodes switchScoreType_(vector<PeptideIdentification>& peptide_ids, const String& id_file_abs_path)
  {
    // Check if score types are valid. TODO
    try
    {
      IDScoreSwitcherAlgorithm switcher;
      Size c = 0;
      switcher.switchToGeneralScoreType(peptide_ids, IDScoreSwitcherAlgorithm::ScoreType::PEP, c);
    }
    catch(Exception::MissingInformation&)
    {
      OPENMS_LOG_FATAL_ERROR << "ProteomicsLFQ expects a Posterior Error Probability score in all Peptide IDs. ID file: " << id_file_abs_path << endl;
      return ExitCodes::INCOMPATIBLE_INPUT_DATA;
    }
    return EXECUTION_OK;
  }

  ExitCodes loadAndCleanupIDFile_(
    const String& id_file_abs_path,
    const String& mz_file,
    const Size& fraction_group,
    const Size& fraction,
    vector<ProteinIdentification>& protein_ids, 
    vector<PeptideIdentification>& peptide_ids,
    set<String>& fixed_modifications,  // adds to
    set<String>& variable_modifications) // adds to
  {
    const String& mz_file_abs_path = File::absolutePath(mz_file);
    IdXMLFile().load(id_file_abs_path, protein_ids, peptide_ids);

    ExitCodes e = checkSingleRunPerID_(protein_ids, id_file_abs_path);
    if (e != EXECUTION_OK) return e;

    e = switchScoreType_(peptide_ids, id_file_abs_path);
    if (e != EXECUTION_OK) return e;
  
    // TODO we could think about removing this limitation 
    IDFilter::keepBestPeptideHits(peptide_ids, false); // strict = false
    IDFilter::removeDecoyHits(peptide_ids);
    IDFilter::removeDecoyHits(protein_ids);
    IDFilter::removeEmptyIdentifications(peptide_ids);
    IDFilter::removeUnreferencedProteins(protein_ids, peptide_ids);

    if (peptide_ids.empty())
    {
      OPENMS_LOG_FATAL_ERROR << "No peptide identifications present after removing decoys " << id_file_abs_path << endl;
      return ExitCodes::INCOMPATIBLE_INPUT_DATA;
    }
 
    // add to the (global) set of fixed and variable modifications
    const vector<String>& var_mods = protein_ids[0].getSearchParameters().variable_modifications;
    const vector<String>& fixed_mods = protein_ids[0].getSearchParameters().fixed_modifications;
    std::copy(var_mods.begin(), var_mods.end(), std::inserter(variable_modifications, variable_modifications.begin())); 
    std::copy(fixed_mods.begin(), fixed_mods.end(), std::inserter(fixed_modifications, fixed_modifications.end())); 

    // delete meta info to free some space
    for (PeptideIdentification & pid : peptide_ids)
    {
      // we currently can't clear the PeptideIdentification meta data
      // because the spectrum_reference is stored in the meta value (which it probably shouldn't)
      // TODO: pid.clearMetaInfo(); if we move it to the PeptideIdentification structure
      for (PeptideHit & ph : pid.getHits())
      {
        // TODO: keep target_decoy information for QC
        ph.clearMetaInfo();
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
      OPENMS_LOG_WARN  << "MS run path not set in ID file: " << id_file_abs_path << endl
                       << "Resetting reference to MS file provided at same input position." << endl;
    }
    else if (id_msfile_ref.size() == 1)
    {
      // Check if the annotated primary MS run filename matches the mzML filename (comparison by base name)
      const String& in_bn = FileHandler::stripExtension(File::basename(mz_file_abs_path));
      const String& id_primaryMSRun_bn = FileHandler::stripExtension(File::basename(id_msfile_ref[0]));

      if (in_bn != id_primaryMSRun_bn)  // mismatch between annotation in ID file and provided mzML file
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
    id_msfile_ref = StringList{mz_file};
    protein_ids[0].setPrimaryMSRunPath(id_msfile_ref);
    protein_ids[0].setMetaValue("fraction_group", fraction_group);
    protein_ids[0].setMetaValue("fraction", fraction);

    // update identifiers to make them unique
    // fixes some bugs related to users splitting the original mzML and id files before running the analysis
    // in that case these files might have the same identifier
    const String old_identifier = protein_ids[0].getIdentifier();
    const String new_identifier = old_identifier + "_" + String(fraction_group) + "F" + String(fraction);
    protein_ids[0].setIdentifier(new_identifier);
    for (PeptideIdentification & p : peptide_ids)
    {
      if (p.getIdentifier() == old_identifier)
      {
        p.setIdentifier(new_identifier);
      }
      else
      {
        OPENMS_LOG_WARN << "Peptide ID identifier found not present in the protein ID" << endl;
      }
    }

    bool missing_spec_ref(false);
    for (const PeptideIdentification & pid : peptide_ids)
    {
      if (!pid.metaValueExists("spectrum_reference") 
        || pid.getMetaValue("spectrum_reference").toString().empty()) 
      {          
        missing_spec_ref = true;
        break;
      }
    }
    // reannotate spectrum references if missing
    if (missing_spec_ref)
    {
      OPENMS_LOG_WARN << "Warning: The identification files don't contain a meta value with the spectrum native id.\n"
                         "OpenMS will try to reannotate them by matching retention times between id and spectra." << endl;

      SpectrumMetaDataLookup::addMissingSpectrumReferences(
        peptide_ids, 
        mz_file_abs_path,
        true);
    }
    return EXECUTION_OK;
  }
 
  ExitCodes quantifyFraction_(
    const pair<unsigned int, std::vector<String> > & ms_files, 
    const map<String, String>& mzfile2idfile, 
    double median_fwhm,
    const multimap<Size, PeptideIdentification> & transfered_ids,
    ConsensusMap & consensus_fraction,
    vector<TransformationDescription> & transformations,
    double& max_alignment_diff,
    set<String>& fixed_modifications,
    set<String>& variable_modifications)
  {
    vector<FeatureMap> feature_maps;
    const Size fraction = ms_files.first;

    const bool is_already_aligned = !transformations.empty();

    // debug output
    writeDebug_("Processing fraction number: " + String(fraction) + "\nFiles: ",  1);
    for (String const & mz_file : ms_files.second) { writeDebug_(mz_file,  1); }

    // for sanity checks we collect the primary MS run basenames as well as the ones stored in the ID files (below)
    StringList id_MS_run_ref;
    StringList in_MS_run = ms_files.second;

    // for each MS file of current fraction (e.g., all MS files that measured the n-th fraction) 
    Size fraction_group{1};
    for (String const & mz_file : ms_files.second)
    { 
      // centroid spectra (if in profile mode) and correct precursor masses
      MSExperiment ms_centroided;    

      {
        ExitCodes e = centroidAndCorrectPrecursors_(mz_file, ms_centroided);
        if (e != EXECUTION_OK) { return e; }
      }

      // load and clean identification data associated with MS run
      vector<ProteinIdentification> protein_ids;
      vector<PeptideIdentification> peptide_ids;
      const String& mz_file_abs_path = File::absolutePath(mz_file);
      const String& id_file_abs_path = File::absolutePath(mzfile2idfile.at(mz_file_abs_path));

      {
        ExitCodes e = loadAndCleanupIDFile_(id_file_abs_path, mz_file, fraction_group, fraction, protein_ids, peptide_ids, fixed_modifications, variable_modifications);
        if (e != EXECUTION_OK) return e;
      }

      StringList id_msfile_ref;
      protein_ids[0].getPrimaryMSRunPath(id_msfile_ref);
      id_MS_run_ref.push_back(id_msfile_ref[0]);
     
      //-------------------------------------------------------------
      // Internal Calibration of spectra peaks and precursor peaks with high-confidence IDs
      //-------------------------------------------------------------
      if (getStringOption_("mass_recalibration") == "true")
      {
        recalibrateMasses_(ms_centroided, peptide_ids, id_file_abs_path);
      }

      vector<ProteinIdentification> ext_protein_ids;
      vector<PeptideIdentification> ext_peptide_ids;

      //////////////////////////////////////////////////////
      // Transfer aligned IDs
      //////////////////////////////////////////////////////
      if (!transfered_ids.empty())
      {
        OPENMS_PRECONDITION(is_already_aligned, "Data has not been aligned.")

        // transform observed IDs and spectra
        MapAlignmentTransformer::transformRetentionTimes(peptide_ids, transformations[fraction_group - 1]);
        MapAlignmentTransformer::transformRetentionTimes(ms_centroided, transformations[fraction_group - 1]);

        // copy the (already) aligned, consensus feature derived ids that are to be transferred to this map to peptide_ids
        auto range = transfered_ids.equal_range(fraction_group - 1);
        for (auto& it = range.first; it != range.second; ++it)
        {
           PeptideIdentification trans = it->second;
           trans.setIdentifier(protein_ids[0].getIdentifier());
           peptide_ids.push_back(trans);
        }
      }

      //////////////////////////////////////////
      // Chromatographic parameter estimation
      //////////////////////////////////////////
      median_fwhm = estimateMedianChromatographicFWHM_(ms_centroided);

      //-------------------------------------------------------------
      // Feature detection
      //-------------------------------------------------------------   
      ///////////////////////////////////////////////

      // Run MTD before FFM

      // create empty feature map and annotate MS file
      FeatureMap seeds;

      StringList sl;
      sl.push_back(mz_file);
      seeds.setPrimaryMSRunPath(sl);

      if (getStringOption_("targeted_only") == "false")
      {
        calculateSeeds_(ms_centroided, seeds, median_fwhm);
        if (debug_level_ > 666)
        {
          FeatureXMLFile().store("debug_seeds_fraction_" + String(ms_files.first) + "_" + String(fraction_group) + ".featureXML", seeds);
        }
      }

      /////////////////////////////////////////////////
      // Run FeatureFinderIdentification

      FeatureMap fm;
      StringList feature_msfile_ref;
      feature_msfile_ref.push_back(mz_file);
      fm.setPrimaryMSRunPath(feature_msfile_ref);

      FeatureFinderIdentificationAlgorithm ffi;
      ffi.getMSData().swap(ms_centroided);
      ffi.getProgressLogger().setLogType(log_type_);

      Param ffi_param = getParam_().copy("PeptideQuantification:", true);
      ffi_param.setValue("detect:peak_width", 5.0 * median_fwhm);
      ffi.setParameters(ffi_param);
      writeDebug_("Parameters passed to FeatureFinderIdentification algorithm", ffi_param, 3);

      FeatureMap tmp = fm;
      ffi.run(peptide_ids, 
        protein_ids, 
        ext_peptide_ids, 
        ext_protein_ids, 
        tmp,
        seeds);

      // TODO: consider moving this to FFid
      // free parts of feature map not needed for further processing (e.g., subfeatures...)
      for (auto & f : tmp)
      {
        //TODO keep FWHM meta value for QC
        f.clearMetaInfo();
        f.setSubordinates({});
        f.setConvexHulls({});
      }

      IDConflictResolverAlgorithm::resolve(tmp, false); // keep only best peptide per feature

      feature_maps.push_back(tmp);
      
      if (debug_level_ > 666)
      {
        FeatureXMLFile().store("debug_fraction_" + String(ms_files.first) + "_" + String(fraction_group) + ".featureXML", feature_maps.back());
      }

      ++fraction_group;
    }

    // Check for common mistake that order of input files have been switched.
    // This is the case if basenames are identical but the order does not match.
    if (!File::validateMatchingFileNames(in_MS_run, id_MS_run_ref, true, true, false)) // only basenames, without extension, only order
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__,
        OPENMS_PRETTY_FUNCTION, "MS run path reference in ID files and spectra filenames match but order differs.");
    }

    //-------------------------------------------------------------
    // Align all features of this fraction (if not already aligned)
    //-------------------------------------------------------------
    if (!is_already_aligned)
    {
      max_alignment_diff = alignAndLink_(
        feature_maps, 
        consensus_fraction, 
        transformations,
        median_fwhm);
    }
    else // Data already aligned. Link with previously determined alignment difference
    {
      link_(feature_maps,
        median_fwhm,
        max_alignment_diff,
        consensus_fraction);
    }

    // add dataprocessing
    if (feature_maps.size() > 1)
    {
      addDataProcessing_(consensus_fraction,
        getProcessingInfo_(DataProcessing::ALIGNMENT));
      addDataProcessing_(consensus_fraction,
        getProcessingInfo_(DataProcessing::FEATURE_GROUPING));
    }

    ////////////////////////////////////////////////////////////
    // Annotate experimental design in consensus map
    ////////////////////////////////////////////////////////////
    Size j(0);
    // for each MS file (as provided in the experimental design)
    for (String const & mz_file : ms_files.second) 
    {
      const Size curr_fraction_group = j + 1;
      consensus_fraction.getColumnHeaders()[j].label = "label-free";
      consensus_fraction.getColumnHeaders()[j].filename = mz_file;
      consensus_fraction.getColumnHeaders()[j].unique_id = feature_maps[j].getUniqueId();
      consensus_fraction.getColumnHeaders()[j].setMetaValue("fraction", fraction);
      consensus_fraction.getColumnHeaders()[j].setMetaValue("fraction_group", curr_fraction_group);
      ++j;
    }

    // assign unique ids
    consensus_fraction.applyMemberFunction(&UniqueIdInterface::setUniqueId);

    // sort list of peptide identifications in each consensus feature by map index
    consensus_fraction.sortPeptideIdentificationsByMapIndex();

    if (debug_level_ >= 666)
    {
      ConsensusXMLFile().store("debug_fraction_" + String(ms_files.first) +  ".consensusXML", consensus_fraction);
      writeDebug_("to produce a consensus map with: " + String(consensus_fraction.getColumnHeaders().size()) + " columns.", 1);
    }

    //-------------------------------------------------------------
    // ID conflict resolution
    //-------------------------------------------------------------
    IDConflictResolverAlgorithm::resolve(consensus_fraction, true);

    //-------------------------------------------------------------
    // ConsensusMap normalization (basic)
    //-------------------------------------------------------------
    if (getStringOption_("out_msstats").empty())  // only normalize if no MSstats output is generated
    {
      ConsensusMapNormalizerAlgorithmMedian::normalizeMaps(
        consensus_fraction, 
        ConsensusMapNormalizerAlgorithmMedian::NM_SCALE, 
        "", 
        "");
    }

    // max_alignment_diff returned by reference
    return EXECUTION_OK;
  }


  ExitCodes inferProteinGroups_(const StringList in_ids,
    const String& in_db,
    const map<String, String>& idfile2mzfile,
    const set<String>& fixed_modifications, 
    vector<ProteinIdentification>& inferred_protein_ids, 
    vector<PeptideIdentification>& inferred_peptide_ids)
  {
    // load the IDs again and merge
    IDMergerAlgorithm merger{String("all_merged")};

    IdXMLFile f;
    
    for (const auto& idfile : in_ids)
    {
      vector<ProteinIdentification> protein_ids;
      vector<PeptideIdentification> peptide_ids;
      f.load(idfile, protein_ids, peptide_ids);

      // Check if score types are valid.
      //TODO do that in the inference algorithms? Epifany only does it for consensusXML
      try
      {
        IDScoreSwitcherAlgorithm switcher;
        Size c = 0;
        switcher.switchToGeneralScoreType(peptide_ids, IDScoreSwitcherAlgorithm::ScoreType::PEP, c);
      }
      catch(Exception::MissingInformation&)
      {
        OPENMS_LOG_FATAL_ERROR <<
          "ProteomicsLFQ expects a Posterior Error Probability score in all Peptide IDs. ID file: "
          << idfile << endl;
        return ExitCodes::INCOMPATIBLE_INPUT_DATA;
      }

      //TODO we could think about removing this limitation
      IDFilter::keepBestPeptideHits(peptide_ids, false); // strict = false

      // reannotate MS run if not present
      StringList id_msfile_ref;
      protein_ids[0].getPrimaryMSRunPath(id_msfile_ref);
      if (id_msfile_ref.empty())
      {
        id_msfile_ref.push_back(idfile2mzfile.at(idfile));
        protein_ids[0].setPrimaryMSRunPath(id_msfile_ref);
      }     
 
      // TODO: Filter for a PSM FDR? Better on an experiment-level though
      merger.insertRuns(std::move(protein_ids), std::move(peptide_ids));
    }

    // For now, we merge all into one. Inference per condition would be another option
    merger.returnResultsAndClear(inferred_protein_ids[0], inferred_peptide_ids);

    if (debug_level_ >= 666)
    {
      IdXMLFile().store("debug_mergedIDs.idXML", inferred_protein_ids, inferred_peptide_ids);
    }

    // since we don't require an index as input but need to calculate e.g., coverage we reindex here (fast)
    if (!in_db.empty())
    {
      PeptideIndexing indexer;
      Param param_pi = indexer.getParameters();
      param_pi.setValue("missing_decoy_action", "silent");
      param_pi.setValue("write_protein_sequence", "true");
      param_pi.setValue("write_protein_description", "true");
      indexer.setParameters(param_pi);

      // stream data in fasta file
      FASTAContainer<TFI_File> fasta_db(in_db);
      PeptideIndexing::ExitCodes indexer_exit = indexer.run(fasta_db, inferred_protein_ids, inferred_peptide_ids);

      if ((indexer_exit != PeptideIndexing::EXECUTION_OK) &&
          (indexer_exit != PeptideIndexing::PEPTIDE_IDS_EMPTY))
      {
        if (indexer_exit == PeptideIndexing::DATABASE_EMPTY)
        {
          return INPUT_FILE_EMPTY;
        }
        else if (indexer_exit == PeptideIndexing::UNEXPECTED_RESULT)
        {
          return UNEXPECTED_RESULT;
        }
        else
        {
          return UNKNOWN_ERROR;
        }
      }
    }

    //-------------------------------------------------------------
    // Protein inference
    //-------------------------------------------------------------
    // TODO: Think about ProteinInference on IDs only merged per condition
    bool groups = getStringOption_("protein_quantification") != "strictly_unique_peptides";
    bool bayesian = getStringOption_("protein_inference") == "bayesian";
    if (!bayesian) // simple aggregation
    {
      BasicProteinInferenceAlgorithm bpia;
      bpia.run(inferred_peptide_ids, inferred_protein_ids);

      if (groups)
      {
        IDBoostGraph ibg{inferred_protein_ids[0], inferred_peptide_ids, 0, false, false};
        ibg.computeConnectedComponents();
        ibg.calculateAndAnnotateIndistProteins(true);
        auto & ipg = inferred_protein_ids[0].getIndistinguishableProteins();
        std::sort(std::begin(ipg), std::end(ipg));
      }
    }
    else // if (bayesian)
    {
      if (!groups) //TODO @julianus: easy to fix. Remove that limitation by adding a bool param.
      {
        throw OpenMS::Exception::InvalidParameter(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "Inference = Bayes currently automatically groups proteins and does not allow for"
          " strictly_unique_peptides during quantification.");
      }
      // Important Note: BayesianProteinInference by default keeps only the best
      // PSM per peptide!
      // TODO maybe allow otherwise!
      BayesianProteinInferenceAlgorithm bayes;
      //bayesian inference automatically annotates groups
      bayes.inferPosteriorProbabilities(inferred_protein_ids, inferred_peptide_ids);
    }

    // if no or only partial grouping was performed, add rest of proteins as singleton groups
    // (assumed by greedy resolution for now)
    inferred_protein_ids[0].fillIndistinguishableGroupsWithSingletons();

    if (debug_level_ >= 666)
    {
      IdXMLFile().store("debug_mergedIDs_inference.idXML", inferred_protein_ids, inferred_peptide_ids);
    }

    // TODO think about order of the next three steps (greedy resolution, FDR calc and filtering)

    // Optional greedy group resolution
    // TODO finish greedy resolution on the new graph structure, since we built it for inference already anyway
    bool greedy_group_resolution = getStringOption_("protein_quantification") == "shared_peptides";
    if (greedy_group_resolution)
    {
      PeptideProteinResolution ppr{};
      ppr.buildGraph(inferred_protein_ids[0], inferred_peptide_ids);
      ppr.resolveGraph(inferred_protein_ids[0], inferred_peptide_ids);
      // TODO maybe move the removal as an option into the PPResolution class
      // TODO add an option to calculate FDR including those "second best protein hits"?
      IDFilter::removeUnreferencedProteins(inferred_protein_ids, inferred_peptide_ids);
      IDFilter::updateProteinGroups(inferred_protein_ids[0].getIndistinguishableProteins(), inferred_protein_ids[0].getHits());
      IDFilter::updateProteinGroups(inferred_protein_ids[0].getProteinGroups(), inferred_protein_ids[0].getHits());
      if (debug_level_ >= 666)
      {
        IdXMLFile().store("debug_mergedIDsGreedyResolved.idXML", inferred_protein_ids, inferred_peptide_ids);
      }
    }

    //-------------------------------------------------------------
    // Protein (and additional peptide?) FDR
    //-------------------------------------------------------------
    const double max_fdr = getDoubleOption_("proteinFDR");
    // Note: actually, when Bayesian inference was performed, only one (best) PSM
    // is left per peptide, so the calculated PSM FDR is equal to a Peptide FDR
    const double max_psm_fdr = getDoubleOption_("psmFDR");
    FalseDiscoveryRate fdr;
    fdr.applyBasic(inferred_protein_ids[0]);
    if (max_psm_fdr < 1.)
    {
      fdr.applyBasic(inferred_peptide_ids);
    }

    if (debug_level_ >= 666)
    {
      // This is needed because we throw out decoy proteins during FDR
      IDFilter::updateProteinReferences(inferred_peptide_ids, inferred_protein_ids, true);
      IDFilter::removeUnreferencedProteins(inferred_protein_ids, inferred_peptide_ids); // if we dont filter peptides for now, we dont need this
      IDFilter::updateProteinGroups(inferred_protein_ids[0].getIndistinguishableProteins(), inferred_protein_ids[0].getHits());
      IDFilter::updateProteinGroups(inferred_protein_ids[0].getProteinGroups(), inferred_protein_ids[0].getHits());

      IdXMLFile().store("debug_mergedIDsGreedyResolvedFDR.idXML", inferred_protein_ids, inferred_peptide_ids);
    }

    // FDR filtering
    if (max_psm_fdr < 1.)
    {
      IDFilter::filterHitsByScore(inferred_peptide_ids, max_psm_fdr);
    }
    if (max_fdr < 1.)
    {
      IDFilter::filterHitsByScore(inferred_protein_ids, max_fdr);
      IDFilter::updateProteinReferences(inferred_peptide_ids, inferred_protein_ids, true);
    }
    if (max_psm_fdr < 1.)
    {
      IDFilter::removeUnreferencedProteins(inferred_protein_ids,
                                           inferred_peptide_ids);
    }
    if (max_fdr < 1. || max_psm_fdr < 1.)
    {
      IDFilter::updateProteinGroups(inferred_protein_ids[0].getIndistinguishableProteins(), inferred_protein_ids[0].getHits());
      IDFilter::updateProteinGroups(inferred_protein_ids[0].getProteinGroups(), inferred_protein_ids[0].getHits());
    }

    if (inferred_protein_ids[0].getHits().empty())
    {
      throw Exception::MissingInformation(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "No proteins left after FDR filtering. Please check the log and adjust your settings.");
    }

    if (debug_level_ >= 666)
    {
      IdXMLFile().store("debug_mergedIDsGreedyResolvedFDRFiltered.idXML", inferred_protein_ids, inferred_peptide_ids);
    }

    // ensure that only one final inference result is generated for now
    assert(inferred_protein_ids.size() == 1);

    // do we only want to keep strictly unique peptides (e.g., no groups)?
    if (!greedy_group_resolution && !groups)
    {
      IDFilter::keepUniquePeptidesPerProtein(inferred_peptide_ids);
      if (debug_level_ >= 666)
      {
        IdXMLFile().store("debug_mergedIDsFDRFilteredStrictlyUniqueResolved.idXML", inferred_protein_ids, inferred_peptide_ids);
      }
    }

    // compute coverage (sequence was annotated during PeptideIndexing)
    inferred_protein_ids[0].computeCoverage(inferred_peptide_ids);

    // TODO: this might not be correct if only the best peptidoform is kept
    // determine observed modifications (exclude fixed mods)
    inferred_protein_ids[0].computeModifications(inferred_peptide_ids, StringList(fixed_modifications.begin(), fixed_modifications.end()));

    return EXECUTION_OK;
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // Parameter handling
    //-------------------------------------------------------------    

    // Read tool parameters
    StringList in = getStringList_("in");
    String out = getStringOption_("out");
    String out_msstats = getStringOption_("out_msstats");
    StringList in_ids = getStringList_("ids");
    String design_file = getStringOption_("design");
    String in_db = getStringOption_("fasta");

    // Validate parameters
    if (in.size() != in_ids.size())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, 
        OPENMS_PRETTY_FUNCTION, "Number of spectra file (" + String(in.size()) + ") must match number of ID files (" + String(in_ids.size()) + ").");
    }

    if (getStringOption_("quantification_method") == "spectral_counting" && !out_msstats.empty())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, 
        OPENMS_PRETTY_FUNCTION, "MSstats export for spectral counting data not supported. Please remove output file.");
    }

    //-------------------------------------------------------------
    // Experimental design: read or generate default
    //-------------------------------------------------------------      
    ExperimentalDesign design;
    if (!design_file.empty())
    { // load from file
      design = ExperimentalDesignFile::load(design_file, false);
      // some sanity checks
      if (design.getNumberOfLabels() != 1)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, 
          OPENMS_PRETTY_FUNCTION, "Experimental design is not label-free as it contains multiple labels.");          
      }
      if (!design.sameNrOfMSFilesPerFraction())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, 
          OPENMS_PRETTY_FUNCTION, "Different number of fractions for different samples provided. This is currently not supported by ProteomicsLFQ.");          
      }
      
      // extract basenames from experimental design and input files
      const auto& pl2fg = design.getPathLabelToFractionGroupMapping(true);      
      set<String> ed_basenames;
      for (const auto& p : pl2fg)
      {
        const String& filename = p.first.first;
        ed_basenames.insert(filename);
      }

      set<String> in_basenames;
      for (Size i = 0; i != in.size(); ++i)
      {
        const String& in_bn = File::basename(in[i]);
        in_basenames.insert(in_bn);
      }

      if (ed_basenames != in_basenames)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__,
          OPENMS_PRETTY_FUNCTION, "Spectra files provided as input need to match the ones in the experimental design file.");
      }
    }
    else
    {
      OPENMS_LOG_INFO << "No experimental design file provided.\n"
                      << "Assuming a label-free experiment without fractionation.\n"
                      << endl;

      // default to unfractionated design
      ExperimentalDesign::MSFileSection msfs;
      Size count{1};
      for (String & s : in)
      {
        ExperimentalDesign::MSFileSectionEntry e;
        e.fraction = 1;
        e.fraction_group = count;
        e.label = 1;
        e.path = s;
        e.sample = count;
        msfs.push_back(e);
      }      
      design.setMSFileSection(msfs);
    }
    std::map<unsigned int, std::vector<String> > frac2ms = design.getFractionToMSFilesMapping();

    for (auto & f : frac2ms)
    {
      writeDebug_("Fraction " + String(f.first) + ":", 10);
      for (const String & s : f.second)
      {
        writeDebug_("MS file: " + s, 10);
      }
    }

    // Map between mzML file and corresponding id file
    // Here we currently assume that these are provided in the exact same order.
    map<String, String> mzfile2idfile = mapMzML2Ids_(in, in_ids);
    map<String, String> idfile2mzfile = mapId2MzMLs_(mzfile2idfile);

    // check if mzMLs in experimental design match to mzMLs passed as in parameter
    for (auto const & ms_files : frac2ms) // for each fraction->ms file(s)
    {      
      for (String const & mz_file : ms_files.second)
      { 
        const String& mz_file_abs_path = File::absolutePath(mz_file);
        if (mzfile2idfile.find(mz_file_abs_path) == mzfile2idfile.end())
        {
          OPENMS_LOG_FATAL_ERROR << "MzML file in experimental design file '"
            << mz_file_abs_path << "'not passed as 'in' parameter.\n" 
            << "Note: relative paths in the experimental design file "
            << "are resolved relative to the design file path. \n"
            << "Use absolute paths or make sure the design file is in "
            << "the same path as the mzML files."
            << endl;
          throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, mz_file_abs_path);
        }
      }
    }

    Param pep_param = getParam_().copy("Posterior Error Probability:", true);
    writeDebug_("Parameters passed to PEP algorithm", pep_param, 3);

    // TODO: inference parameter

    Param pq_param = getParam_().copy("ProteinQuantification:", true);
    writeDebug_("Parameters passed to PeptideAndProteinQuant algorithm", pq_param, 3);


    Param com_param = getParam_().copy("algorithm:common:", true);
    writeDebug_("Common parameters passed to both sub-algorithms (mtd and epd)", com_param, 3);

    set<String> fixed_modifications, variable_modifications;

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
      double median_fwhm(0);
      for (auto const & ms_files : frac2ms) // for each fraction->ms file(s)
      {      
        ConsensusMap consensus_fraction; // quantitative result for this fraction identifier
        vector<TransformationDescription> transformations; // filled by RT alignment
        double max_alignment_diff(0.0);

        ExitCodes e = quantifyFraction_(
          ms_files, 
          mzfile2idfile,
          median_fwhm, 
          multimap<Size, PeptideIdentification>(),
          consensus_fraction, 
          transformations,  // transformations are empty, will be filled by alignment
          max_alignment_diff,  // max_alignment_diff not yet determined, will be filled by alignment
          fixed_modifications, 
          variable_modifications);

        if (e != EXECUTION_OK) { return e; }
        
        if (getStringOption_("transfer_ids") != "false")
        {  
          OPENMS_LOG_INFO << "Transferring identification data between runs of the same fraction." << endl;
          // needs to occur in >= 50% of all runs for transfer
          const Size min_occurrance = (ms_files.second.size() + 1) / 2;
          multimap<Size, PeptideIdentification> transfered_ids = transferIDsBetweenSameFraction_(consensus_fraction, min_occurrance);
          consensus_fraction.clear();

          // The transferred IDs were calculated on the aligned data
          // So we make sure we use the aligned IDs and peak maps in the re-quantification step
          e = quantifyFraction_(
            ms_files, 
            mzfile2idfile, 
            median_fwhm, 
            transfered_ids, 
            consensus_fraction, 
            transformations,  // transformations as determined by alignment
            max_alignment_diff, // max_alignment_error as determined by alignment
            fixed_modifications,
            variable_modifications);

          OPENMS_POSTCONDITION(!consensus_fraction.empty(), "ConsensusMap of fraction empty after ID transfer.!");
          if (e != EXECUTION_OK) { return e; }
        }
        consensus.appendColumns(consensus_fraction);  // append consensus map calculated for this fraction number
      }  // end of scope of fraction related data

      consensus.sortByPosition();
      consensus.sortPeptideIdentificationsByMapIndex();

      if (debug_level_ >= 666)
      {
        ConsensusXMLFile().store("debug_after_normalization.consensusXML", consensus);
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
      for (auto const & ms_files : frac2ms) // for each fraction->ms file(s)
      {
        const Size& fraction = ms_files.first;

        // debug output
        writeDebug_("Processing fraction number: " + String(fraction) + "\nFiles: ",  1);
        for (String const & mz_file : ms_files.second) { writeDebug_(mz_file,  1); }

        // for sanity checks we collect the primary MS run basenames as well as the ones stored in the ID files (below)
        StringList id_MS_run_ref;
        StringList in_MS_run = ms_files.second;

        // for each MS file of current fraction (e.g., all MS files that measured the n-th fraction) 
        Size fraction_group{1};
        for (String const & mz_file : ms_files.second)
        { 
          // load and clean identification data associated with MS run
          vector<ProteinIdentification> protein_ids;
          vector<PeptideIdentification> peptide_ids;
          const String& mz_file_abs_path = File::absolutePath(mz_file);
          const String& id_file_abs_path = File::absolutePath(mzfile2idfile.at(mz_file_abs_path));

          {
            ExitCodes e = loadAndCleanupIDFile_(id_file_abs_path, mz_file, fraction_group, fraction, protein_ids, peptide_ids, fixed_modifications, variable_modifications);
            if (e != EXECUTION_OK) return e;
          }

          StringList id_msfile_ref;
          protein_ids[0].getPrimaryMSRunPath(id_msfile_ref);
          id_MS_run_ref.push_back(id_msfile_ref[0]);

          // append to consensus map
          all_protein_ids.emplace_back(std::move(protein_ids[0]));
          all_peptide_ids.insert(all_peptide_ids.end(), 
            std::make_move_iterator(peptide_ids.begin()), 
            std::make_move_iterator(peptide_ids.end()));
        }

        ////////////////////////////////////////////////////////////
        // Annotate experimental design in consensus map
        ////////////////////////////////////////////////////////////
        Size j(0);
        // for each MS file (as provided in the experimental design)
        for (String const & mz_file : ms_files.second) 
        {
          const Size curr_fraction_group = j + 1;
          consensus.getColumnHeaders()[run_index].label = "label-free";
          consensus.getColumnHeaders()[run_index].filename = mz_file;
          consensus.getColumnHeaders()[run_index].unique_id = 1 + run_index;
          consensus.getColumnHeaders()[run_index].setMetaValue("fraction", fraction);
          consensus.getColumnHeaders()[run_index].setMetaValue("fraction_group", curr_fraction_group);
          ++j;
          ++run_index;
        }
      }
    }

    //-------------------------------------------------------------
    // ID related algorithms
    // TODO we could switch to work on the IDs in ConsensusXML
    //  but not all algorithms are available on Cons.Maps yet.
    //-------------------------------------------------------------

    vector<ProteinIdentification> inferred_protein_ids{1};
    vector<PeptideIdentification> inferred_peptide_ids;
    ExitCodes e = inferProteinGroups_(in_ids, in_db, idfile2mzfile, fixed_modifications, inferred_protein_ids, inferred_peptide_ids);
    if (e != EXECUTION_OK) return e;
   
    // references from PSM to Protein that got removed during inference need to be also removed in the consensus map
    {
      set<String> protein_inferred;
      // determine all inferred proteins
      for (auto& p : inferred_protein_ids)
      {
        for (auto& ph : p.getHits())
        {
          protein_inferred.insert(ph.getAccession());
        }
      }

      // update assigned and unassigned peptide identifications
      for (auto& c : consensus)
      {
        auto& pids = c.getPeptideIdentifications();
        for (auto& pid : pids)
        {
          auto& hits = pid.getHits();
          for (auto& ph : hits)
          {
            std::vector<PeptideEvidence> pes = ph.getPeptideEvidences();
            pes.erase(std::remove_if(pes.begin(), 
                                pes.end(),
                                [&protein_inferred](PeptideEvidence& x){ return protein_inferred.find(x.getProteinAccession()) == protein_inferred.end(); }),
                 pes.end());
            ph.setPeptideEvidences(std::move(pes));
          }
        // erase hits without reference to proteins
        hits.erase(std::remove_if(hits.begin(), hits.end(),
                              [](PeptideHit& x){ return x.getPeptideEvidences().empty(); }),
               hits.end());
        }
        // erase peptide ids without peptide hits
        pids.erase(std::remove_if(pids.begin(), 
                         pids.end(),
                         [](PeptideIdentification& x){ return x.getHits().empty(); }),
               pids.end());
      }

      auto& pids = consensus.getUnassignedPeptideIdentifications();
      for (auto& pid : pids)
      {
        auto& hits = pid.getHits();
        for (auto& ph : hits)
        {
          std::vector<PeptideEvidence> pes = ph.getPeptideEvidences();
          pes.erase(std::remove_if(pes.begin(), 
                              pes.end(),
                              [&protein_inferred](PeptideEvidence& x){ return protein_inferred.find(x.getProteinAccession()) == protein_inferred.end(); }),
               pes.end());
          ph.setPeptideEvidences(std::move(pes));
        }

        // erase hits without reference to proteins
        hits.erase(std::remove_if(hits.begin(), hits.end(),
                              [](PeptideHit& x){ return x.getPeptideEvidences().empty(); }),
               hits.end());
      }

     // erase peptide ids without peptide hits
     pids.erase(std::remove_if(pids.begin(), 
                         pids.end(),
                         [](PeptideIdentification& x){ return x.getHits().empty(); }),
          pids.end());
    }

    // clean up references (assigned and unassigned)
    IDFilter::removeUnreferencedProteins(consensus, true);

    //-------------------------------------------------------------
    // Peptide quantification
    //-------------------------------------------------------------
    PeptideAndProteinQuant quantifier;

    if (getStringOption_("quantification_method") == "feature_intensity")
    {
      quantifier.setParameters(pq_param);
      quantifier.readQuantData(consensus, design);
    }
    else if (getStringOption_("quantification_method") == "spectral_counting")
    {
      pq_param.setValue("average", "sum"); 
      pq_param.setValue("top", 0); // all 
      pq_param.setValue("consensus:normalize", "false");
      quantifier.setParameters(pq_param);

      quantifier.readQuantData(
       consensus.getProteinIdentifications(),
       consensus.getUnassignedPeptideIdentifications(),
       design);
    }

    quantifier.quantifyPeptides(inferred_peptide_ids);

    //-------------------------------------------------------------
    // Protein quantification
    //-------------------------------------------------------------

    // Should always be there by now, even if just singletons (TODO a bit of a waste then, though)
    if (inferred_protein_ids[0].getIndistinguishableProteins().empty())
    {
      throw Exception::MissingInformation(
       __FILE__, 
       __LINE__, 
       OPENMS_PRETTY_FUNCTION, 
       "No information on indistinguishable protein groups found.");
    }

    quantifier.quantifyProteins(inferred_protein_ids[0]);
    auto const & protein_quants = quantifier.getProteinResults();
    if (protein_quants.empty())
    {        
     OPENMS_LOG_WARN << "Warning: No proteins were quantified." << endl;
    }

    if (debug_level_ >= 666)
    {
      IdXMLFile().store("debug_quant.idXML", inferred_protein_ids, inferred_peptide_ids);
    }
    //-------------------------------------------------------------
    // Export of MzTab file as final output
    //-------------------------------------------------------------

    // Annotate quants to protein(groups) for easier export in mzTab
    // Note: we keep protein groups that have not been quantified
    PeptideAndProteinQuant::annotateQuantificationsToProteins(protein_quants, inferred_protein_ids[0], design.getNumberOfFractionGroups(), false);

    if (debug_level_ >= 666)
    {
      IdXMLFile().store("debug_quant_annotated.idXML", inferred_protein_ids, inferred_peptide_ids);
    }

    // insert inference information as first protein identification
    auto& proteins = consensus.getProteinIdentifications();
    proteins.insert(proteins.begin(), inferred_protein_ids[0]);

    // For correctness we would need to set the run reference in the pepIDs of the consensusXML all to the first run then
    // And probably make sure that peptides that correspond to filtered out proteins are not producing errors
    // e.g. by removing them with a Filter beforehand.

    consensus.resolveUniqueIdConflicts(); // TODO: find out if this is still needed
    if (!getStringOption_("out_cxml").empty())
    {
      // Note: idXML and consensusXML doesn't support writing quantification at protein groups
      // (they are neverless stored and passed to mzTab for proper export)
      ConsensusXMLFile().store(getStringOption_("out_cxml"), consensus);
    }

    // Fill MzTab with meta data and quants annotated in identification data structure
    const bool report_unmapped(true);
    const bool report_unidentified_features(false);

    MzTab m = MzTab::exportConsensusMapToMzTab(
      consensus, 
      String("null"),
      true,
      report_unidentified_features, 
      report_unmapped,
      "Export from ProteomicsLFQ workflow in OpenMS.");
    MzTabFile().store(out, m);

    if (!out_msstats.empty())
    {
      MSstatsFile msstats;
      // TODO: add a helper method to quickly check if experimental design file contain the right columns (and put this at start of tool)

      // shrink protein runs to the one containing the inference data
      consensus.getProteinIdentifications().resize(1);

      msstats.storeLFQ(
        out_msstats, 
        consensus, 
        design, 
        StringList(), 
        false, 
        "MSstats_BioReplicate", 
        "MSstats_Condition", 
        "max");
    }

    return EXECUTION_OK;
  }
};

int main(int argc, const char ** argv)
{
  ProteomicsLFQ tool;
  return tool.main(argc, argv);
}

/// @endcond
