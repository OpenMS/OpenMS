// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>

#include <OpenMS/METADATA/MSQuantifications.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzQuantMLFile.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringCentroided.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringProfile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexClustering.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexSatelliteCentroided.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>
#include <OpenMS/COMPARISON/CLUSTERING/GridBasedCluster.h>
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
using namespace boost::math;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_FeatureFinderMultiplex FeatureFinderMultiplex
  @brief Detects peptide pairs in LC-MS data and determines their relative abundance.
<CENTER>
  <table>
    <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FeatureFinderMultiplex \f$ \longrightarrow \f$</td>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
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
    TOPPBase("FeatureFinderMultiplex", "Determination of peak ratios in LC-MS data", true)
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
    FeatureXMLFile file;
    file.store(filename, map);
  }
  
  /**
   * @brief Write consensus map to consensusXML file.
   *
   * @param filename    name of consensusXML file
   * @param map    consensus map for output
   */
  void writeConsensusMap_(const String& filename, ConsensusMap& map) const
  {     
    ConsensusXMLFile file;
    for (auto & ch : map.getColumnHeaders())
    {
      ch.second.filename = getStringOption_("in");
    }
    file.store(filename, map);
  }
  
  /**
   * @brief Write blacklist to mzML file.
   *
   * @param filename    name of mzML file
   * @param blacklist    blacklist for output
   */
  void writeBlacklist_(const String& filename, const MSExperiment& blacklist) const
  {    
    MzMLFile file;
    file.store(filename, blacklist);
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
    MzMLFile file;
    MSExperiment exp;
    
    // only read MS1 spectra
    std::vector<int> levels;
    levels.push_back(1);
    file.getOptions().setMSLevels(levels);
    
    LOG_DEBUG << "Loading input..." << endl;
    file.setLogType(log_type_);
    file.load(in_, exp);

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
      writeFeatureMap_(out_, algorithm.getFeatureMap());
    }
    if (!(out_multiplets_.empty()))
    {
      writeConsensusMap_(out_multiplets_, algorithm.getConsensusMap());
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

//@endcond
