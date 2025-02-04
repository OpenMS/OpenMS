// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------


#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>

#include <OpenMS/FORMAT/FileHandler.h>

#include <typeinfo>

using namespace OpenMS;
using namespace std;

/**
@page TOPP_SpectraFilterThresholdMower SpectraFilterThresholdMower

@brief Removes all peaks below an intensity threshold.

<CENTER>
<table>
<tr>
<th ALIGN = "center"> pot. predecessor tools </td>
<td VALIGN="middle" ROWSPAN=2> &rarr; SpectraFilter &rarr;</td>
<th ALIGN = "center"> pot. successor tools </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating on MS peak data @n (in mzML format)</td>
</tr>
</table>
</CENTER>

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_SpectraFilterThresholdMower.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_SpectraFilterThresholdMower.html
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSpectraFilterThresholdMower :
  public TOPPBase
{
public:
  TOPPSpectraFilterThresholdMower() :
    TOPPBase("SpectraFilterThresholdMower", "Applies thresholdfilter to peak spectra.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output file ");
    setValidFormats_("out", ListUtils::create<String>("mzML"));

    // register one section for each algorithm
    registerSubsection_("algorithm", "Algorithm parameter subsection.");

  }

  Param getSubsectionDefaults_(const String & /*section*/) const override
  {
    return ThresholdMower().getParameters();
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    //input/output files
    String in(getStringOption_("in"));
    String out(getStringOption_("out"));

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    PeakMap exp;
    FileHandler().loadExperiment(in, exp, {FileTypes::MZML}, log_type_);

    //-------------------------------------------------------------
    // if meta data arrays are present, remove them and warn
    //-------------------------------------------------------------
    if (exp.clearMetaDataArrays())
    {
      writeLogWarn_("Warning: Spectrum meta data arrays cannot be sorted. They are deleted.");
    }

    //-------------------------------------------------------------
    // filter
    //-------------------------------------------------------------
    Param filter_param = getParam_().copy("algorithm:", true);
    writeDebug_("Used filter parameters", filter_param, 3);

    ThresholdMower filter;
    filter.setParameters(filter_param);
    filter.filterPeakMap(exp);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    //annotate output with data processing info
    addDataProcessing_(exp, getProcessingInfo_(DataProcessing::FILTERING));

    FileHandler().storeExperiment(out, exp, {FileTypes::MZML}, log_type_);

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPSpectraFilterThresholdMower tool;
  return tool.main(argc, argv);
}

/// @endcond

