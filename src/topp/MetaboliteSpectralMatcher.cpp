// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/FORMAT/MSPGenericFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_MetaboliteSpectralMatcher MetaboliteSpectralMatcher

@brief MetaboliteSpectralMatcher identifies small molecules from tandem MS spectra using a spectral library.

<CENTER>
<table>
<tr>
<th ALIGN = "center"> pot. predecessor tools </td>
<td VALIGN="middle" ROWSPAN=3> &rarr; MetaboliteSpectralMatcher &rarr;</td>
<th ALIGN = "center"> pot. successor tools </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> processing in R </td>
</tr>
</table>
</CENTER>

By default, MS2 spectra with similar precursor mass are merged before comparison with database spectra, for example when a mass at the beginning of the peak and on the peak apex is selected twice as precursor.
Merging can also have disadvantages, for example, for isobaric or isomeric compounds that have similar/same masses but can have different retention times and MS2 spectra.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_MetaboliteSpectralMatcher.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_MetaboliteSpectralMatcher.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMetaboliteSpectralMatcher :
        public TOPPBase
{
public:
  TOPPMetaboliteSpectralMatcher() :
      TOPPBase("MetaboliteSpectralMatcher", "Perform a spectral library search.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input spectra.");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerInputFile_("database", "<file>", "", "Default spectral database.", true);
    setValidFormats_("database", {"mzML", "msp", "mgf"});
    registerOutputFile_("out", "<file>", "", "mzTab file");
    setValidFormats_("out", ListUtils::create<String>("mzTab"));
    registerOutputFile_("out_spectra", "<file>", "", "Output spectra as mzML file. Can be useful to inspect the peak map after spectra merging.", false);
    setValidFormats_("out_spectra", ListUtils::create<String>("mzML"));

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
    String out_spectra = getStringOption_("out_spectra");

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    FileHandler mz_file;
    std::vector<Int> ms_level = {2};
    mz_file.getOptions().setMSLevels(ms_level);

    PeakMap ms_peakmap;
    mz_file.loadExperiment(in, ms_peakmap, {FileTypes::MZML});

    if (ms_peakmap.empty())
    {
      OPENMS_LOG_WARN << "The input file does not contain any MS2/fragment spectra.";
      return INCOMPATIBLE_INPUT_DATA;
    }

    MzTab mztab_output;
    MzTabFile mztab_outfile;

    //-------------------------------------------------------------
    // get parameters
    //-------------------------------------------------------------

    Param msm_param = getParam_().copy("algorithm:", true);
    writeDebug_("Parameters passed to MetaboliteSpectralMatcher", msm_param, 3);

    //-------------------------------------------------------------
    // load database
    //-------------------------------------------------------------
    PeakMap spec_db;
    FileHandler().loadExperiment(spec_db_filename, spec_db, {FileTypes::MSP, FileTypes::MZML, FileTypes::MGF});

    if (spec_db.empty())
    {
      OPENMS_LOG_WARN << "The spectral library does not contain any spectra.";
      return INCOMPATIBLE_INPUT_DATA;
    }

    //-------------------------------------------------------------
    // run spectral library search
    //-------------------------------------------------------------
    MetaboliteSpectralMatching msm;
    msm.setParameters(msm_param);
    msm.run(ms_peakmap, spec_db, mztab_output, out_spectra);

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

