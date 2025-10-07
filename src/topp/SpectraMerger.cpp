// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Andreas Bertsch, Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/PROCESSING/SPECTRAMERGING/SpectraMerger.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_SpectraMerger SpectraMerger

@brief Allows to add up several spectra.

<center>
<table>
<tr>
<th ALIGN = "center"> pot. predecessor tools </td>
<td VALIGN="middle" ROWSPAN=2> &rarr; SpectraMerger &rarr;</td>
<th ALIGN = "center"> pot. successor tools </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating on MS peak data @n (in mzML format) </td>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating on MS peak data @n (in mzML format)</td>
</tr>
</table>
</center>

@experimental This TOPP-tool is not well tested and not all features might be properly implemented and tested!

This tool can add several consecutive scans, increasing S/N ratio (for MS1 and above) or merge scans which stem from similar precursors (for MS2 and above).

In any case, the number of scans will be reduced.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_SpectraMerger.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_SpectraMerger.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSpectraMerger :
  public TOPPBase
{
public:
  TOPPSpectraMerger() :
    TOPPBase("SpectraMerger", "Merges spectra (each MS level separately), increasing S/N ratios.")
  {

  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input mzML file.");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Output mzML file with merged spectra.");
    setValidFormats_("out", ListUtils::create<String>("mzML"));

    registerStringOption_("merging_method", "<method>", "average_gaussian", "Method of merging which should be used.", false);
    setValidStrings_("merging_method", ListUtils::create<String>("average_gaussian,average_tophat,precursor_method,block_method"));

    registerSubsection_("algorithm", "Algorithm section for merging spectra");
  }

  Param getSubsectionDefaults_(const String & /*section*/) const override
  {
    return SpectraMerger().getParameters();
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in(getStringOption_("in"));
    String out(getStringOption_("out"));
    String merging_method(getStringOption_("merging_method"));

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    FileHandler fh;
    FileTypes::Type in_type = fh.getType(in);

    PeakMap exp;
    fh.loadExperiment(in, exp, {in_type}, log_type_);
    exp.sortSpectra();
    exp.updateRanges();

    auto levels = exp.getMSLevels();
    if (levels.empty()) throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, levels.size());
    int min_ms_level = levels.front();
    int max_ms_level = levels.back();
    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    SpectraMerger merger;
    merger.setLogType(log_type_);
    merger.setParameters(getParam_().copy("algorithm:", true));
    if (merging_method == "precursor_method")
    {
      merger.mergeSpectraPrecursors(exp);
    }
    else if (merging_method == "block_method")
    {
      merger.mergeSpectraBlockWise(exp);
    }
    else if (merging_method == "average_gaussian")
    {
      int ms_level = merger.getParameters().getValue("average_gaussian:ms_level");
      if (ms_level == 0)
      {
        for (int tmp_ms_level = min_ms_level; tmp_ms_level <= max_ms_level; tmp_ms_level++)
        {
          merger.average(exp, "gaussian", tmp_ms_level);
        }
      }
      else
      {
        merger.average(exp, "gaussian", ms_level);
      }
    }
    else if (merging_method == "average_tophat")
    {
      int ms_level = merger.getParameters().getValue("average_tophat:ms_level");
      if (ms_level == 0)
      {
        for (int tmp_ms_level = min_ms_level; tmp_ms_level <= max_ms_level; tmp_ms_level++)
        {
          merger.average(exp, "tophat", tmp_ms_level);
        }
      }
      else
      {
        merger.average(exp, "tophat", ms_level);
      }
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------


    fh.storeExperiment(out, exp, {}, log_type_);

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPSpectraMerger tool;
  return tool.main(argc, argv);
}

/// @endcond
