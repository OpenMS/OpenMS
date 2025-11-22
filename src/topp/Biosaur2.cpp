// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2023.
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
// $Authors: Mark Ivanov, Timo Sachsenberg $
// --------------------------------------------------------------------------

/**
 * @page TOPP_Biosaur2 Biosaur2
 *
 * @brief Feature detection for LC-MS1 data
 *
 * This TOPP tool is a C++ reimplementation of the Biosaur2 feature detection algorithm.
 * It detects peptide features in centroided LC-MS1 data by:
 * 1. Grouping peaks across scans into "hills"
 * 2. Splitting hills at valley points
 * 3. Detecting isotope patterns
 * 4. Calculating feature properties (m/z, RT, intensity, charge)
 *
 * Reference:
 * Abdrakhimov, et al. Biosaur: An open-source Python software for liquid
 * chromatography-mass spectrometry peptide feature detection with ion mobility support.
 * https://doi.org/10.1002/rcm.9045
 *
 * <B>The command line parameters of this tool are:</B>
 * @verbinclude TOPP_Biosaur2.cli
 * <B>INI file documentation of this tool:</B>
 * @htmlinclude TOPP_Biosaur2.html
 */

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FEATUREFINDER/Biosaur2Algorithm.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <vector>

using namespace OpenMS;
using namespace std;

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPBiosaur2 final :
  public TOPPBase
{
public:
  TOPPBiosaur2() :
    TOPPBase("Biosaur2", "Feature detection for LC-MS1 data", false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    Biosaur2Algorithm algo_defaults;
    const Param& defaults = algo_defaults.getDefaults();

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

    MSExperiment exp;
    MzMLFile mzml_file;
    mzml_file.load(in, exp);

    Biosaur2Algorithm algorithm;
    Param algo_param = getParam_().copySubset(Biosaur2Algorithm().getDefaults());
    algorithm.setParameters(algo_param);

    FeatureMap feature_map;
    vector<Biosaur2Algorithm::Hill> hills;
    vector<Biosaur2Algorithm::PeptideFeature> peptide_features;
    algorithm.run(exp, feature_map, hills, peptide_features);
    addDataProcessing_(feature_map, getProcessingInfo_(DataProcessing::DATA_PROCESSING));

    FeatureXMLFile feature_file;
    feature_file.store(out, feature_map);
    OPENMS_LOG_INFO << "Wrote " << peptide_features.size() << " features to: " << out << endl;

    if (write_hills_flag || !out_hills.empty())
    {
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
      algorithm.writeHills(hills, hills_file);
    }

    if (!out_tsv.empty())
    {
      algorithm.writeTSV(peptide_features, out_tsv);
    }

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPBiosaur2 tool;
  return tool.main(argc, argv);
}

/// @endcond
