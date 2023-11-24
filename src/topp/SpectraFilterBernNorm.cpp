// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------


#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FILTERING/TRANSFORMERS/BernNorm.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <typeinfo>

using namespace OpenMS;
using namespace std;

/**
  @page TOPP_SpectraFilterBernNorm SpectraFilterBernNorm

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
  @verbinclude TOPP_SpectraFilterBernNorm.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_SpectraFilterBernNorm.html
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSpectraFilterBernNorm :
  public TOPPBase
{
public:
  TOPPSpectraFilterBernNorm() :
    TOPPBase("SpectraFilterBernNorm", "Scales and filters spectra according using the Bern norm.")
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
    return BernNorm().getParameters();
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

    BernNorm filter;
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
  TOPPSpectraFilterBernNorm tool;
  return tool.main(argc, argv);
}

/// @endcond

