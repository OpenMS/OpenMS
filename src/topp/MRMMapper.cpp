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

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/ANALYSIS/TARGETED/MRMMapping.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_MRMMapper MRMMapper

  @brief MRMMapper maps measured chromatograms (mzML) and the transitions used (TraML).

  <CENTER>
      <table>
          <tr>
              <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
              <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ MRMMapper \f$ \longrightarrow \f$</td>
              <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
          </tr>
          <tr>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_FileFilter </td>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathAnalyzer </td>
          </tr>
          <tr>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref UTILS_MRMTransitionGroupPicker </td>
          </tr>
      </table>
  </CENTER>
 
  This tool reads an mzML containing chromatograms (presumably measured on an
  SRM instrument) and a TraML file that contains the data that was used to
  generate the instrument method to measure said data. It then maps the
  transitions in the TraML file to the chromatograms found in the mzML file
  and stores the chromatograms annotated with meta-data from the TraML file.
  Thus, the  output chromatograms are an annotated copy of the input
  chromatograms with native id, precursor information and peptide sequence (if
  available) annotated in the chromatogram files.

  The algorithm tries to match a given set of chromatograms and targeted
  assays. It iterates through all the chromatograms retrieves one or more
  matching targeted assay for the chromatogram. By default, the algorithm
  assumes that a 1:1 mapping exists. If a chromatogram cannot be mapped
  (does not have a corresponding assay) the algorithm issues a warning, the
  user can specify that the program should abort in such a case (see
  error_on_unmapped).
      
  If multiple mapping is enabled (see map_multiple_assays parameter)
  then each mapped assay will get its own chromatogram that contains the
  same raw data but different meta-annotation. This *can* be useful if the
  same transition is used to monitor multiple analytes but may also
  indicate a problem with too wide mapping tolerances.

  The thus mapped mzML file can then be used in a downstream analysis.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_MRMMapper.cli

  <B>The algorithm parameters for the Analyzer filter are:</B>
  @htmlinclude TOPP_MRMMapper.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPMRMMapper 
  : public TOPPBase
{

public:

  TOPPMRMMapper() :
    TOPPBase("MRMMapper", "MRMMapper maps measured chromatograms (mzML) and the transitions used (TraML)", true)
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file containing chromatograms (converted mzXML file)");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("tr", "<file>", "", "transition file");
    setValidFormats_("tr", ListUtils::create<String>("TraML"));

    registerOutputFile_("out", "<file>", "", "Output file containing mapped chromatograms");
    setValidFormats_("out", ListUtils::create<String>("mzML"));

    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String& name) const
  {
    if (name == "algorithm")
    {
      return MRMMapping().getDefaults();
    }
    else
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown subsection", name);
    }
  }

  ExitCodes main_(int, const char **)
  {
    String in = getStringOption_("in");
    String tr_file = getStringOption_("tr");
    String out = getStringOption_("out");

    OpenMS::TargetedExperiment targeted_exp;
    OpenMS::PeakMap chromatogram_map;
    OpenMS::PeakMap output;

    TraMLFile().load(tr_file, targeted_exp);
    MzMLFile().load(in, chromatogram_map);

    Param param = getParam_().copy("algorithm:", true);

    MRMMapping mrmm;
    mrmm.setParameters(param);
    mrmm.mapExperiment(chromatogram_map, targeted_exp, output);

    // add all data processing information to all the chromatograms
    DataProcessing dp_ = getProcessingInfo_(DataProcessing::FORMAT_CONVERSION);
    DataProcessingPtr dp = boost::shared_ptr<DataProcessing>(new DataProcessing(dp_));
    std::vector<MSChromatogram > chromatograms = output.getChromatograms();
    for (Size i=0; i<chromatograms.size(); ++i)
    {
      chromatograms[i].getDataProcessing().push_back(dp);
    }
    output.setChromatograms(chromatograms);

    MzMLFile().store(out, output);
    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{

  TOPPMRMMapper tool;
  return tool.main(argc, argv);
}

/// @endcond
