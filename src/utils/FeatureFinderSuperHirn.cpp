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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/Feature.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSH.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_FeatureFinderSuperHirn FeatureFinderSuperHirn

  A feature finder based on the original SuperHirn codebase.

  Proteomics. 2007 Oct;7(19):3470-80.
  SuperHirn - a novel tool for high resolution LC-MS-based peptide/protein profiling.
  Mueller LN, Rinner O, Schmidt A, Letarte S, Bodenmiller B, Brusniak MY, Vitek O, Aebersold R, Mller M.

  see FeatureFinderAlgorithmSH.h for a detailed description of the algorithm.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

typedef FeatureFinderAlgorithmSH FFSH;

class TOPPFeatureFinderSH :
  public TOPPBase
{
public:
  TOPPFeatureFinderSH() :
    TOPPBase("FeatureFinderSuperHirn", "Finds mass spectrometric features in mass spectra.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input profile data file ");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output peak file ");
    setValidFormats_("out", ListUtils::create<String>("featureXML"));

    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    return FFSH().getDefaults();
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out = getStringOption_("out");

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    MzMLFile mzMLFile;
    mzMLFile.setLogType(log_type_);
    PeakMap input;
    mzMLFile.getOptions().addMSLevel(1);
    mzMLFile.load(in, input);

    if (input.empty())
    {
      LOG_WARN << "The given file does not contain any conventional peak data, but might"
                  " contain chromatograms. This tool currently cannot handle them, sorry.";
      return INCOMPATIBLE_INPUT_DATA;
    }

    //check if spectra are sorted
    for (Size i = 0; i < input.size(); ++i)
    {
      if (!input[i].isSorted())
      {
        writeLog_("Error: Not all spectra are sorted according to peak m/z positions. Use FileFilter to sort the input!");
        return INCOMPATIBLE_INPUT_DATA;
      }
    }

    //-------------------------------------------------------------
    // pick
    //-------------------------------------------------------------
    FeatureMap output;

    FeatureFinder ff;
    Param param = getParam_().copy("algorithm:", true);

    FFSH ffsh;
    ffsh.setParameters(param);
    ffsh.setData(input, output, ff);
    ffsh.run();

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    //annotate output with data processing info
    addDataProcessing_(output, getProcessingInfo_(DataProcessing::PEAK_PICKING));
    addDataProcessing_(output, getProcessingInfo_(DataProcessing::QUANTITATION));
    output.ensureUniqueId();
    for (Size i = 0; i < output.size(); i++)
    {
      output[i].ensureUniqueId();
    }
    StringList ms_runs;
    input.getPrimaryMSRunPath(ms_runs);
    output.setPrimaryMSRunPath(ms_runs);
    FeatureXMLFile().store(out, output);

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPFeatureFinderSH tool;
  return tool.main(argc, argv);
}

/// @endcond
