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
// $Authors: Erhan Kenar $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/MzTabFile.h>

#include <OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
        @page UTILS_MetaboliteSpectralMatcher MetaboliteSpectralMatcher

        @brief MetaboliteSpectralMatcher identify small molecules from tandem MS spectra.

        <CENTER>
        <table>
        <tr>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ MetaboliteSpectralMatcher \f$ \longrightarrow \f$</td>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref PeakPickerHiRes </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> processing in R </td>
        </tr>
        </table>
        </CENTER>

        <B>The command line parameters of this tool are:</B>
        @verbinclude TOPP_MetaboliteSpectralMatcher.cli


        MetaboliteSpectralMatcher matches spectra from a spectral library with tandem MS spectra.

        <B>The command line parameters of this tool are:</B>
        @verbinclude UTILS_MetaboliteSpectralMatcher.cli
        <B>INI file documentation of this tool:</B>
        @htmlinclude UTILS_MetaboliteSpectralMatcher.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMetaboliteSpectralMatcher :
        public TOPPBase
{
public:
  TOPPMetaboliteSpectralMatcher() :
      TOPPBase("MetaboliteSpectralMatcher", "Perform a spectral library search.", false)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input spectra.");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerInputFile_("database", "<file>", "CHEMISTRY/MetaboliteSpectralDB.mzML", "Default spectral database.", false);
    setValidFormats_("database", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "mzTab file");
    setValidFormats_("out", ListUtils::create<String>("tsv"));

    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    return MetaboliteSpectralMatching().getDefaults();
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String database = getStringOption_("database");
    String spec_db_filename(database);

    // default path? retrieve file path in share folder
    if (database == "CHEMISTRY/MetaboliteSpectralDB.mzML")
    {
      // throws Exception::FileNotFound if file does not exist
      spec_db_filename = File::find("CHEMISTRY/MetaboliteSpectralDB.mzML");
    }

    String out = getStringOption_("out");

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    MzMLFile mz_file;
    mz_file.setLogType(log_type_);
    std::vector<Int> ms_level(1,2);
    mz_file.getOptions().setMSLevels(ms_level);

    PeakMap ms_peakmap;
    mz_file.load(in, ms_peakmap);

    if (ms_peakmap.empty())
    {
      LOG_WARN << "The input file does not contain any spectra.";
      return INCOMPATIBLE_INPUT_DATA;
    }

    MzTab mztab_output;
    MzTabFile mztab_outfile;

    //-------------------------------------------------------------
    // get parameters
    //-------------------------------------------------------------

    Param ams_param = getParam_().copy("algorithm:", true);
    writeDebug_("Parameters passed to MetaboliteSpectralMatcher", ams_param, 3);

    //-------------------------------------------------------------
    // load database
    //-------------------------------------------------------------

    PeakMap spec_db;
    mz_file.load(spec_db_filename, spec_db);

    if (spec_db.empty())
    {
      LOG_WARN << "The spectral library does not contain any spectra.";
      return INCOMPATIBLE_INPUT_DATA;
    }

    //-------------------------------------------------------------
    // run spectral library search
    //-------------------------------------------------------------
    MetaboliteSpectralMatching ams;
    ams.setParameters(ams_param);
    ams.run(ms_peakmap, spec_db, mztab_output);

    //-------------------------------------------------------------
    // store results
    //-------------------------------------------------------------
    mztab_outfile.store(out, mztab_output);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPMetaboliteSpectralMatcher tool;
  return tool.main(argc, argv);
}

/// @endcond

