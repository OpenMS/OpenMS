// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
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
          <th ALIGN = "center"> potential predecessor tools </td>
          <td VALIGN="middle" ROWSPAN=3> &rarr; MRMMapper &rarr;</td>
          <th ALIGN = "center"> potential successor tools </td>
      </tr>
      <tr>
          <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_FileFilter </td>
          <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathAnalyzer </td>
      </tr>
      <tr>
          <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MRMTransitionGroupPicker </td>
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

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file containing chromatograms (converted mzXML file)");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("tr", "<file>", "", "transition file");
    setValidFormats_("tr", ListUtils::create<String>("traML"));

    registerOutputFile_("out", "<file>", "", "Output file containing mapped chromatograms");
    setValidFormats_("out", ListUtils::create<String>("mzML"));

    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String& name) const override
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

  ExitCodes main_(int, const char **) override
  {
    String in = getStringOption_("in");
    String tr_file = getStringOption_("tr");
    String out = getStringOption_("out");

    OpenMS::TargetedExperiment targeted_exp;
    OpenMS::PeakMap chromatogram_map;
    OpenMS::PeakMap output;

    FileHandler().loadTransitions(tr_file, targeted_exp, {FileTypes::TRAML});
    FileHandler().loadExperiment(in, chromatogram_map, {FileTypes::MZML});

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

    FileHandler().storeExperiment(out, output, {FileTypes::MZML});
    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{

  TOPPMRMMapper tool;
  return tool.main(argc, argv);
}

/// @endcond
