// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors:  Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzQuantMLFile.h>
#include <OpenMS/METADATA/MSQuantifications.h>
#include <OpenMS/SYSTEM/File.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
 @page TOPP_FeatureFinderCentroided FeatureFinderCentroided

 @brief The feature detection application for quantitation (centroided).

<CENTER>
 <table>
  <tr>
   <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
   <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FeatureFinderCentroided \f$ \longrightarrow \f$</td>
   <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
  </tr>
  <tr>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerWavelet </td>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureLinkerUnlabeled @n (or another feature grouping tool) </td>
  </tr>
  <tr>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_SeedListGenerator </td>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MapAlignerPoseClustering @n (or another alignment tool) </td>
  </tr>
 </table>
</CENTER>

 Reference:\n
 Weisser <em>et al.</em>: <a href="http://dx.doi.org/10.1021/pr300992u">An automated pipeline for high-throughput label-free quantitative proteomics</a> (J. Proteome Res., 2013, PMID: 23391308).

 This module identifies "features" in a LC/MS map. By feature, we understand a peptide in a MS sample that
 reveals a characteristic isotope distribution. The algorithm
 computes positions in rt and m/z dimension and a charge estimate
 of each peptide.

 The algorithm identifies pronounced regions of the data around so-called <tt>seeds</tt>.
 In the next step, we iteratively fit a model of the isotope profile and the retention time to
 these data points. Data points with a low probability under this model are removed from the
 feature region. The intensity of the feature is then given by the sum of the data points included
 in its regions.

 How to find suitable parameters and details of the different algorithms implemented are described
 in the @ref TOPP_example_featuredetection "TOPP tutorial".

 Specialized tools are available for some experimental techniques: @ref TOPP_IsobaricAnalyzer.

 <B>The command line parameters of this tool are:</B>
 @verbinclude TOPP_FeatureFinderCentroided.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_FeatureFinderCentroided.html

 For the parameters of the algorithm section see the algorithms documentation: @n
  @ref OpenMS::FeatureFinderAlgorithmPicked "centroided" @n

 In the following table you can find example values of the most important parameters for
 different instrument types. @n These parameters are not valid for all instruments of that type,
 but can be used as a starting point for finding suitable parameters.

 <b>'centroided' algorithm</b>:
 <table>
  <tr>
   <td>&nbsp;</td>
   <td><b>Q-TOF</b></td>
   <td><b>LTQ Orbitrap</b></td>
  </tr>
  <tr>
   <td><b>intensity:bins</b></td>
   <td>10</td>
   <td>10</td>
  </tr>
  <tr>
   <td><b>mass_trace:mz_tolerance</b></td>
   <td>0.02</td>
   <td>0.004</td>
  </tr>
  <tr>
   <td><b>isotopic_pattern:mz_tolerance</b></td>
   <td>0.04</td>
   <td>0.005</td>
  </tr>
 </table>

 For the @em centroided algorithm centroided data is needed. In order to create centroided data from profile data use the @ref TOPP_PeakPickerWavelet.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderCentroided :
  public TOPPBase
{
public:
  TOPPFeatureFinderCentroided() :
    TOPPBase("FeatureFinderCentroided", "Detects two-dimensional features in LC-MS data.")
  {}

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("featureXML"));
    registerInputFile_("seeds", "<file>", "", "User specified seed list", false);
    setValidFormats_("seeds", ListUtils::create<String>("featureXML"));

    registerOutputFile_("out_mzq", "<file>", "", "Optional output file of MzQuantML.", false, true);
    setValidFormats_("out_mzq", ListUtils::create<String>("mzq"));

    addEmptyLine_();

    registerSubsection_("algorithm", "Algorithm section");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    return FeatureFinder().getParameters(FeatureFinderAlgorithmPicked::getProductName());
  }

  ExitCodes main_(int, const char**) override
  {
    //input file names
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String out_mzq = getStringOption_("out_mzq");

    //prevent loading of fragment spectra
    PeakFileOptions options;
    options.setMSLevels(vector<Int>(1, 1));

    //reading input data
    MzMLFile f;
    f.getOptions() = options;
    f.setLogType(log_type_);

    PeakMap exp;
    f.load(in, exp);
    exp.updateRanges();

    if (exp.getSpectra().empty())
    {
      throw OpenMS::Exception::FileEmpty(__FILE__, __LINE__, __FUNCTION__, "Error: No MS1 spectra in input file.");
    }

    // determine type of spectral data (profile or centroided)
    SpectrumSettings::SpectrumType  spectrum_type = exp[0].getType();

    if (spectrum_type == SpectrumSettings::PROFILE)
    {
      if (!getFlag_("force"))
      {
        throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __FUNCTION__, "Error: Profile data provided but centroided spectra expected. To enforce processing of the data set the -force flag.");
      }
    }

    //load seeds
    FeatureMap seeds;
    if (getStringOption_("seeds") != "")
    {
      FeatureXMLFile().load(getStringOption_("seeds"), seeds);
    }

    //setup of FeatureFinder
    FeatureFinder ff;
    ff.setLogType(log_type_);

    // A map for the resulting features
    FeatureMap features;
    StringList ms_runs;
    exp.getPrimaryMSRunPath(ms_runs);
    features.setPrimaryMSRunPath(ms_runs);

    // get parameters specific for the feature finder
    Param feafi_param = getParam_().copy("algorithm:", true);
    writeDebug_("Parameters passed to FeatureFinder", feafi_param, 3);

    // Apply the feature finder
    ff.run(FeatureFinderAlgorithmPicked::getProductName(), exp, features, feafi_param, seeds);
    features.applyMemberFunction(&UniqueIdInterface::setUniqueId);

    // DEBUG
    if (debug_level_ > 10)
    {
      FeatureMap::Iterator it;
      for (it = features.begin(); it != features.end(); ++it)
      {
        if (!it->isMetaEmpty())
        {
          vector<String> keys;
          it->getKeys(keys);
          LOG_INFO << "Feature " << it->getUniqueId() << endl;
          for (Size i = 0; i < keys.size(); i++)
          {
            LOG_INFO << "  " << keys[i] << " = " << it->getMetaValue(keys[i]) << endl;
          }
        }
      }
    }

    //-------------------------------------------------------------
    // writing files
    //-------------------------------------------------------------

    //annotate output with data processing info
    addDataProcessing_(features, getProcessingInfo_(DataProcessing::QUANTITATION));

    // write features to user specified output file
    FeatureXMLFile map_file;

    // Remove detailed convex hull information and subordinate features
    // (unless requested otherwise) to reduce file size of feature files
    // unless debugging is turned on.
    if (debug_level_ < 5)
    {
      FeatureMap::Iterator it;
      for (it = features.begin(); it != features.end(); ++it)
      {
        it->getConvexHull().expandToBoundingBox();
        for (Size i = 0; i < it->getConvexHulls().size(); ++i)
        {
          it->getConvexHulls()[i].expandToBoundingBox();
        }
        it->getSubordinates().clear();
      }
    }

    map_file.store(out, features);

    if (!out_mzq.trim().empty())
    {
      std::vector<DataProcessing> tmp;
      for (Size i = 0; i < exp[0].getDataProcessing().size(); i++)
      {
        tmp.push_back(*exp[0].getDataProcessing()[i].get());
      }
      MSQuantifications msq(features, exp.getExperimentalSettings(), tmp );
      msq.assignUIDs();
      MzQuantMLFile file;
      file.store(out_mzq, msq);
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPFeatureFinderCentroided tool;
  return tool.main(argc, argv);
}

/// @endcond
