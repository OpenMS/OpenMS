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
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmIsotopeWavelet.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_FeatureFinderIsotopeWavelet FeatureFinderIsotopeWavelet

    @brief The feature detection application for quantitation.

<CENTER>
    <table>
        <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB" ROWSPAN=1> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FeatureFinderIsotopeWavelet \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_NoiseFilterSGolay </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MapAlignerPoseClustering @n (or another alignment tool) </td>
        </tr>
        <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_NoiseFilterGaussian </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureLinkerUnlabeled @n (or another feature grouping tool) </td>
        </tr>
    </table>
</CENTER>

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

    @note that the wavelet transform is very slow on high-resolution spectra (i.e. FT, Orbitrap). We recommend
    to use a noise or intensity filter to remove spurious points first and to speed-up the feature detection process.

    Specialized tools are available for some experimental techniques: @ref TOPP_IsobaricAnalyzer.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_FeatureFinderIsotopeWavelet.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_FeatureFinderIsotopeWavelet.html

    For the parameters of the algorithm section see the algorithms documentation: @n
    @ref OpenMS::FeatureFinderAlgorithmIsotopeWavelet "isotope_wavelet" @n

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderIsotopeWavelet :
  public TOPPBase
{
public:
  TOPPFeatureFinderIsotopeWavelet() :
    TOPPBase("FeatureFinderIsotopeWavelet", "Detects two-dimensional features in LC-MS data.")
  {}

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("featureXML"));

    addEmptyLine_();
    registerSubsection_("algorithm", "Algorithm section");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    return FeatureFinder().getParameters(FeatureFinderAlgorithmIsotopeWavelet::getProductName());
  }

  ExitCodes main_(int, const char**) override
  {
    //input and output file names ..
    String in = getStringOption_("in");
    String out = getStringOption_("out");

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

    //no seeds supported
    FeatureMap seeds;

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
    ff.run(FeatureFinderAlgorithmIsotopeWavelet::getProductName(), exp, features, feafi_param, seeds);
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
    map_file.store(out, features);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPFeatureFinderIsotopeWavelet tool;
  return tool.main(argc, argv);
}

/// @endcond
