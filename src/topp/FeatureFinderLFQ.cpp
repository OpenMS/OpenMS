// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Mark Ivanov, Timo Sachsenberg $
// --------------------------------------------------------------------------

/**
 * @page TOPP_FeatureFinderLFQ FeatureFinderLFQ
 *
 * @brief Feature detection for LC-MS1 data with ion mobility support (EXPERIMENTAL)
 *
 * <B>Note: This tool is experimental and under active development. The interface and behavior may change.</B>
 *
 * This TOPP tool is a C++ reimplementation of the Biosaur2 feature detection algorithm.
 * It detects peptide features in centroided LC-MS1 data (with optional profile mode support) through:
 * 1. Grouping peaks across scans into "hills" (continuous m/z traces)
 * 2. Splitting hills at valley points to separate co-eluting species
 * 3. Detecting isotope patterns based on expected mass differences and intensity correlations
 * 4. Calculating comprehensive feature properties (m/z, RT, intensity, charge state)
 *
 * <B>Key Features:</B>
 * - FAIMS compensation voltage grouping for FAIMS-enabled instruments
 * - Ion mobility-aware processing for PASEF/TIMS data (2D centroiding in m/z and ion mobility space)
 * - TOF-specific intensity filtering for time-of-flight instruments
 * - Automatic mass calibration to improve detection accuracy
 * - Profile mode support via PeakPickerHiRes centroiding
 * - Export to featureXML and Biosaur2-compatible TSV formats
 *
 * The tool closely mirrors the Python reference implementation to ensure reproducible results
 * and exposes all core parameters through the INI file for fine-tuning. Besides the mandatory
 * featureXML output, optional TSV exports for both the peptide features and raw hills can be
 * enabled for quality control and downstream analysis.
 *
 * <B>Reference:</B>
 * Abdrakhimov, et al. Biosaur: An open-source Python software for liquid chromatography-mass
 * spectrometry peptide feature detection with ion mobility support.
 * Rapid Communications in Mass Spectrometry, 2022. https://doi.org/10.1002/rcm.9045
 *
 * <B>The command line parameters of this tool are:</B>
 * @verbinclude TOPP_FeatureFinderLFQ.cli
 * <B>INI file documentation of this tool:</B>
 * @htmlinclude TOPP_FeatureFinderLFQ.html
 */

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FEATUREFINDER/Biosaur2Algorithm.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/SYSTEM/File.h>
#include <vector>

using namespace OpenMS;
using namespace std;

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

/**
 * @brief TOPP wrapper that exposes the Biosaur2 feature detection pipeline on the command line.
 *
 * The tool handles mzML I/O, forwards the full set of @ref OpenMS::Biosaur2Algorithm parameters and optionally writes TSV exports for
 * features and hills. All algorithmic heavy lifting is performed by the reusable Biosaur2Algorithm class so the TOPP wrapper remains
 * focused on parameter handling and reporting.
 */
class TOPPFeatureFinderLFQ final :
  public TOPPBase
{
public:
  TOPPFeatureFinderLFQ() :
    TOPPBase("FeatureFinderLFQ", "Feature detection for LC-MS1 data (EXPERIMENTAL)", false)
  {
  }

protected:
  /// Declare command-line options and forward the algorithm defaults into the INI namespace.
  void registerOptionsAndFlags_() override
  {
    const Param& defaults = algorithm_.getDefaults();

    registerInputFile_("in", "<file>", "", "Input mzML file (centroided data)");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerOutputFile_("out", "<file>", "", "Output featureXML file");
    setValidFormats_("out", ListUtils::create<String>("featureXML"));

    registerOutputFile_("out_tsv", "<file>", "", "Optional: output TSV file (Biosaur2 format)", false);
    setValidFormats_("out_tsv", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_hills", "<file>", "", "Optional: write detected hills to TSV", false);
    setValidFormats_("out_hills", ListUtils::create<String>("tsv"));

    registerFlag_("write_hills", "Force writing of hills file even if no output path was provided", false);

    registerFullParam_(defaults);
  }

  /**
   * @brief Execute the Biosaur2 workflow for the configured inputs/outputs.
   *
   * The function loads the requested mzML file, instantiates Biosaur2Algorithm with the user-defined parameters and takes care of
   * exporting featureXML, TSV and hill diagnostics. Meta information such as primary MS run paths and DataProcessing entries are annotated.
   */
  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String out_tsv = getStringOption_("out_tsv");
    String out_hills = getStringOption_("out_hills");
    bool write_hills_flag = getFlag_("write_hills");

    Param algo_param = getParam_().copySubset(algorithm_.getDefaults());
    algorithm_.setParameters(algo_param);

    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    StopWatch stopwatch;

    progresslogger.startProgress(0, 1, "Loading input mzML");
    stopwatch.start();
    MzMLFile mzml_file;
    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(1); // only load MS1 level for feature finding
    mzml_file.setOptions(options);
    mzml_file.load(in, algorithm_.getMSData());
    progresslogger.setProgress(1);
    progresslogger.endProgress();
    stopwatch.stop();
    OPENMS_LOG_INFO << "Loaded input file in " << stopwatch.toString() << endl;

    FeatureMap feature_map;
    vector<Biosaur2Algorithm::Hill> hills;
    vector<Biosaur2Algorithm::PeptideFeature> peptide_features;
    stopwatch.reset();
    progresslogger.startProgress(0, 1, "Preprocessing and feature finding");
    stopwatch.start();
    algorithm_.run(feature_map, hills, peptide_features);
    progresslogger.setProgress(1);
    progresslogger.endProgress();
    stopwatch.stop();
    String primary_path = getFlag_("test") ? ("file://" + File::basename(in)) : in;
    feature_map.setPrimaryMSRunPath({primary_path}, algorithm_.getMSData());
    addDataProcessing_(feature_map, getProcessingInfo_(DataProcessing::QUANTITATION));
    OPENMS_LOG_INFO << "Preprocessing and feature finding took " << stopwatch.toString() << endl;

    stopwatch.reset();
    progresslogger.startProgress(0, 1, "Writing featureXML output");
    stopwatch.start();
    FeatureXMLFile feature_file;
    feature_file.store(out, feature_map);
    progresslogger.setProgress(1);
    progresslogger.endProgress();
    stopwatch.stop();
    OPENMS_LOG_INFO << "Wrote " << peptide_features.size() << " features to: " << out << endl;
    OPENMS_LOG_INFO << "Writing featureXML took " << stopwatch.toString() << endl;

    if (write_hills_flag || !out_hills.empty())
    {
      stopwatch.reset();
      String hills_file = out_hills;
      if (hills_file.empty())
      {
        String base = out;
        Size dot_pos = base.find_last_of('.');
        if (dot_pos != String::npos)
        {
          base = base.substr(0, dot_pos);
        }
        hills_file = base + ".hills.tsv";
      }
      progresslogger.startProgress(0, 1, "Writing hills TSV");
      stopwatch.start();
      algorithm_.writeHills(hills, hills_file);
      progresslogger.setProgress(1);
      progresslogger.endProgress();
      stopwatch.stop();
      OPENMS_LOG_INFO << "Writing hills TSV took " << stopwatch.toString() << endl;
    }

    if (!out_tsv.empty())
    {
      stopwatch.reset();
      progresslogger.startProgress(0, 1, "Writing feature TSV");
      stopwatch.start();
      algorithm_.writeTSV(peptide_features, out_tsv);
      progresslogger.setProgress(1);
      progresslogger.endProgress();
      stopwatch.stop();
      OPENMS_LOG_INFO << "Writing feature TSV took " << stopwatch.toString() << endl;
    }

    return EXECUTION_OK;
  }

private:
  Biosaur2Algorithm algorithm_;
};

int main(int argc, const char** argv)
{
  TOPPFeatureFinderLFQ tool;
  return tool.main(argc, argv);
}

/// @endcond
