// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse, Timo Sachsenberg, Samuel Wein, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/ML/REGRESSION/LinearRegression.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>

#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/ML/REGRESSION/LinearRegression.h>

#include <OpenMS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>
#include <OpenMS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteringCentroided.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteringProfile.h>
#include <OpenMS/FEATUREFINDER/MultiplexClustering.h>
#include <OpenMS/FEATUREFINDER/MultiplexSatelliteCentroided.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>
#include <OpenMS/ML/CLUSTERING/GridBasedCluster.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>

//Contrib includes
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <QDir>

//std includes
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <limits>
#include <locale>
#include <iomanip>

using namespace std;
using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_FeatureFinderMultiplex FeatureFinderMultiplex
@brief Detects peptide pairs in LC-MS data and determines their relative abundance.
<CENTER>
  <table>
    <tr>
      <th ALIGN = "center"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=3> &rarr; FeatureFinderMultiplex &rarr;</td>
      <th ALIGN = "center"> pot. successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileConverter </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_IDMapper</td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileFilter </td>
    </tr>
  </table>
</CENTER>
FeatureFinderMultiplex is a tool for the fully automated analysis of quantitative proteomics data. It detects pairs of isotopic envelopes with fixed m/z separation. It requires no prior sequence identification of the peptides. In what follows we outline the algorithm.
<b>Algorithm</b>
The algorithm is divided into three parts: filtering, clustering and linear fitting, see Fig. (d), (e) and (f). In the following discussion let us consider a particular mass spectrum at retention time 1350 s, see Fig. (a). It contains a peptide of mass 1492 Da and its 6 Da heavier labelled counterpart. Both are doubly charged in this instance. Their isotopic envelopes therefore appear at 746 and 749 in the spectrum. The isotopic peaks within each envelope are separated by 0.5. The spectrum was recorded at finite intervals. In order to read accurate intensities at arbitrary m/z we spline-fit over the data, see Fig. (b).
We would like to search for such peptide pairs in our LC-MS data set. As a warm-up let us consider a standard intensity cut-off filter, see Fig. (c). Scanning through the entire m/z range (red dot) only data points with intensities above a certain threshold pass the filter. Unlike such a local filter, the filter used in our algorithm takes intensities at a range of m/z positions into account, see Fig. (d). A data point (red dot) passes if
- all six intensities at m/z, m/z+0.5, m/z+1, m/z+3, m/z+3.5 and m/z+4 lie above a certain threshold,
- the intensity profiles in neighbourhoods around all six m/z positions show a good correlation and
- the relative intensity ratios within a peptide agree up to a factor with the ratios of a theoretic averagine model.
Let us now filter not only a single spectrum but all spectra in our data set. Data points that pass the filter form clusters in the t-m/z plane, see Fig. (e). Each cluster corresponds to the mono-isotopic mass trace of the lightest peptide of a SILAC pattern. We now use hierarchical clustering methods to assign each data point to a specific cluster. The optimum number of clusters is determined by maximizing the silhouette width of the partitioning. Each data point in a cluster corresponds to three pairs of intensities (at [m/z, m/z+3], [m/z+0.5, m/z+3.5] and [m/z+1, m/z+4]). A plot of all intensity pairs in a cluster shows a clear linear correlation, see Fig. (f). Using linear regression we can determine the relative amounts of labelled and unlabelled peptides in the sample.
@image html SILACAnalyzer_algorithm.png
<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_FeatureFinderMultiplex.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_FeatureFinderMultiplex.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderMultiplex :
  public TOPPBase
{
private:

  // input and output files
  String in_;
  String out_;
  String out_multiplets_;
  String out_blacklist_;

public:
  TOPPFeatureFinderMultiplex() :
    TOPPBase("FeatureFinderMultiplex", "Determination of peak ratios in LC-MS data")
  {
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "LC-MS dataset in either centroid or profile mode");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Output file containing the individual peptide features.", false);
    setValidFormats_("out", ListUtils::create<String>("featureXML"));
    registerOutputFile_("out_multiplets", "<file>", "", "Optional output file containing all detected peptide groups (i.e. peptide pairs or triplets or singlets or ..). The m/z-RT positions correspond to the lightest peptide in each group.", false, true);
    setValidFormats_("out_multiplets", ListUtils::create<String>("consensusXML"));
    registerOutputFile_("out_blacklist", "<file>", "", "Optional output file containing all peaks which have been associated with a peptide feature (and subsequently blacklisted).", false, true);
    setValidFormats_("out_blacklist", ListUtils::create<String>("mzML"));

    registerFullParam_(FeatureFinderMultiplexAlgorithm().getDefaults());
  }

  /**
   * @brief process parameters of 'input/output' section
   */
  void getParameters_in_out_()
  {
    in_ = getStringOption_("in");
    out_ = getStringOption_("out");
    out_multiplets_ = getStringOption_("out_multiplets");
    out_blacklist_ = getStringOption_("out_blacklist");
  }

  /**
   * @brief Write feature map to featureXML file.
   *
   * @param filename    name of featureXML file
   * @param map    feature map for output
   */
  void writeFeatureMap_(const String& filename, FeatureMap& map) const
  {
    FileHandler().storeFeatures(filename, map, {FileTypes::FEATUREXML});
  }

  /**
   * @brief Write consensus map to consensusXML file.
   *
   * @param filename    name of consensusXML file
   * @param map    consensus map for output
   */
  void writeConsensusMap_(const String& filename, ConsensusMap& map) const
  {
    for (auto & ch : map.getColumnHeaders())
    {
      ch.second.filename = getStringOption_("in");
    }
    FileHandler().storeConsensusFeatures(filename, map, {FileTypes::CONSENSUSXML});
  }

  /**
   * @brief Write blacklist to mzML file.
   *
   * @param filename    name of mzML file
   * @param blacklist    blacklist for output
   */
  void writeBlacklist_(const String& filename, const MSExperiment& blacklist) const
  {
    FileHandler().storeExperiment(filename, blacklist, {FileTypes::MZML});
  }

  /**
   * @brief determine the number of samples
   * for example n=2 for SILAC, or n=1 for simple feature detection
   *
   * @param labels    string describing the labels
   */
  static size_t numberOfSamples(String labels)
  {
    // samples can be deliminated by any kind of brackets
    labels.substitute("(", "[");
    labels.substitute("{", "[");
    size_t n = std::count(labels.begin(), labels.end(), '[');
    // if no labels are specified, we have n=1 for simple feature detection
    if (n == 0)
    {
      n = 1;
    }

    return n;
  }


  ExitCodes main_(int, const char**) override
  {

    /**
     * handle parameters
     */
    getParameters_in_out_();

    if ((out_.empty()) && (out_multiplets_.empty()))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Strings for all output files are empty. Please specify at least one output file.");
    }

    /**
     * load input
     */
    FileHandler file;
    MSExperiment exp;

    // only read MS1 spectra
    std::vector<int> levels;
    levels.push_back(1);
    file.getOptions().setMSLevels(levels);

    OPENMS_LOG_DEBUG << "Loading input..." << endl;
    file.loadExperiment(in_, exp, {FileTypes::MZML}, log_type_);

    FeatureFinderMultiplexAlgorithm algorithm;
    // pass only relevant parameters to the algorithm and set the log type
    Param params = getParam_();
    params.remove("in");
    params.remove("out");
    params.remove("out_multiplets");
    params.remove("out_blacklist");
    params.remove("log");
    params.remove("debug");
    params.remove("threads");
    params.remove("no_progress");
    params.remove("force");
    params.remove("test");
    algorithm.setParameters(params);
    algorithm.setLogType(this->log_type_);
    // run feature detection algorithm
    algorithm.run(exp, true);

    // write feature map, consensus maps and blacklist
    if (!(out_.empty()))
    {
      FeatureMap& feature_map = algorithm.getFeatureMap();
      feature_map.setPrimaryMSRunPath({in_}, exp);
      writeFeatureMap_(out_, feature_map);
    }
    if (!(out_multiplets_.empty()))
    {
      ConsensusMap& consensus_map = algorithm.getConsensusMap();
      StringList ms_run_paths(numberOfSamples(params.getValue("algorithm:labels").toString()), in_);
      consensus_map.setPrimaryMSRunPath(ms_run_paths, exp);
      writeConsensusMap_(out_multiplets_, consensus_map);
    }
    if (!(out_blacklist_.empty()))
    {
      writeBlacklist_(out_blacklist_, algorithm.getBlacklist());
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPFeatureFinderMultiplex tool;
  return tool.main(argc, argv);
}

///@endcond
