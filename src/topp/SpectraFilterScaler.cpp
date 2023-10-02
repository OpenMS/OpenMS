// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------


#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FILTERING/TRANSFORMERS/Scaler.h>

#include <OpenMS/FORMAT/MzMLFile.h>

#include <typeinfo>

using namespace OpenMS;
using namespace std;

/**
  @page TOPP_SpectraFilterScaler SpectraFilterScaler

  @brief Filters the top Peaks in the given spectra according to a given schema/thresholdset

  <CENTER>
  <table>
  <tr>
  <th ALIGN = "center"> pot. predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=2> &rarr; SpectraFilter &rarr;</td>
  <th ALIGN = "center"> pot. successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerWavelet </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating on MS peak data @n (in mzML format)</td>
  </tr>
  </table>
  </CENTER>

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_SpectraFilterScaler.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_SpectraFilterScaler.html
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSpectraFilterScaler :
  public TOPPBase
{
public:
  TOPPSpectraFilterScaler() :
    TOPPBase("SpectraFilterScaler", "Applies thresholdfilter to peak spectra.")
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
    return Scaler().getParameters();
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
    MzMLFile f;
    f.setLogType(log_type_);
    f.load(in, exp);

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

    Scaler filter;
    filter.setParameters(filter_param);
    filter.filterPeakMap(exp);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    //annotate output with data processing info
    addDataProcessing_(exp, getProcessingInfo_(DataProcessing::FILTERING));

    f.store(out, exp);

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPSpectraFilterScaler tool;
  return tool.main(argc, argv);
}

/// @endcond

