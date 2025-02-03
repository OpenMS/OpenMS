// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_DeMeanderize DeMeanderize

@brief Repairs MALDI experiments which were spotted line by line.



<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_DeMeanderize.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_DeMeanderize.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDeMeanderize :
  public TOPPBase
{
public:
  TOPPDeMeanderize() :
    TOPPBase("DeMeanderize", "Orders the spectra of MALDI spotting plates correctly.")
  {

  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<mzML-file>", "", "Input experiment file, containing the wrongly sorted spectra.");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<mzML-file>", "", "Output experiment file with correctly sorted spectra.");
    setValidFormats_("out", ListUtils::create<String>("mzML"));
    registerIntOption_("num_spots_per_row", "<integer>", 48, "Number of spots in one column, until next row is spotted.", false);
    setMinInt_("num_spots_per_row", 1);
    registerDoubleOption_("RT_distance", "<integer>", 1.0, "RT distance between two spots which is used to calculated pseudo RT.", false, true);
    setMinFloat_("RT_distance", 0.0);
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in(getStringOption_("in"));
    String out(getStringOption_("out"));
    Size num_spots_per_row(getIntOption_("num_spots_per_row"));
    double RT_distance(getDoubleOption_("RT_distance"));

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    PeakMap exp;
    FileHandler().loadExperiment(in, exp, {FileTypes::MZML}, log_type_);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    ProgressLogger pl;
    pl.setLogType(log_type_);
    pl.startProgress(0, exp.size(), "Assigning pseudo RTs.");
    Size num_ms1(0), num_ms1_base(0), row_counter(0);
    bool row_to_reverse(false);
    double actual_RT(0);
    for (Size i = 0; i != exp.size(); ++i)
    {
      pl.setProgress(i);
      if (row_to_reverse)
      {
        actual_RT = (double)(num_ms1_base + (num_spots_per_row - row_counter)) * RT_distance;
        writeDebug_("RT=" + String(actual_RT) + " (modified, row_counter=" + String(row_counter) + ")", 1);
      }
      else
      {
        actual_RT = (double)num_ms1 * RT_distance;
        writeDebug_("RT=" + String(actual_RT), 1);
      }

      exp[i].setRT(actual_RT);

      if (exp[i].getMSLevel() == 1)
      {
        if (++row_counter >= num_spots_per_row)
        {
          row_counter = 0;
          if (row_to_reverse)
          {
            row_to_reverse = false;
          }
          else
          {
            row_to_reverse = true;
          }
        }
        ++num_ms1;
        if (!row_to_reverse)
        {
          num_ms1_base = num_ms1;
        }
      }
    }
    pl.endProgress();

    // sort the spectra according to their new RT
    exp.sortSpectra();

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------


    FileHandler().storeExperiment(out, exp, {FileTypes::MZML}, log_type_);

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPDeMeanderize tool;
  return tool.main(argc, argv);
}

/// @endcond
