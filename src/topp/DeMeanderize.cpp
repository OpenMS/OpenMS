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

@brief Repairs MALDI experiments which were spotted line by line in a meandering pattern.

<CENTER>
<table>
  <tr>
    <th ALIGN = "center"> pot pot pot </th>
  </tr>
  <tr>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1>
    @image html DeMeanderize.png
    </td>
  </tr>
</table>
</CENTER>

<B>Problem Description:</B>

MALDI (Matrix-Assisted Laser Desorption/Ionization) spotting robots often apply samples to
target plates in a "meandering" or "snake-like" pattern to optimize spotting efficiency. This means:
- Row 1 is spotted from left to right: spots 1, 2, 3, ..., N
- Row 2 is spotted from right to left: spots N+1, N+2, ..., 2N (but in reverse physical order!)
- Row 3 is spotted from left to right again: spots 2N+1, 2N+2, ..., 3N
- And so on...

When the mass spectrometer reads these spots in physical order (left to right for all rows),
the spectra from alternating rows have incorrect retention time assignments, appearing in 
reverse chronological order compared to when they were spotted.

<B>Solution:</B>

DeMeanderize corrects this by:
1. Identifying MS1 spectra that belong to reversed rows
2. Recalculating their pseudo-retention times (RT) to reflect the actual spotting order
3. Sorting all spectra by their corrected RT values

<B>Algorithm Details:</B>

The tool processes spectra sequentially and assigns pseudo-RT values:
- For normal rows (first, third, fifth, etc.):
  @code
  RT = spectrum_number * RT_distance
  @endcode

- For reversed rows (second, fourth, sixth, etc.):
  @code
  RT = (row_start_number + (num_spots_per_row - position_in_row)) * RT_distance
  @endcode

The algorithm maintains:
- A counter tracking total MS1 spectra processed
- The position within the current row (0 to num_spots_per_row-1)
- A flag indicating whether the current row should be reversed

After processing @p num_spots_per_row MS1 spectra, the algorithm switches to the next row
and toggles the reversal flag.

<B>Parameters:</B>
- @p num_spots_per_row: Number of spots in each row of the plate (default: 48)
- @p RT_distance: The pseudo-RT spacing between adjacent spots (default: 1.0)

<B>Example:</B>

For a plate with 4 spots per row and RT_distance=1.0:
@code
Original physical order:  1  2  3  4  8  7  6  5  9 10 11 12
Spotting order:           1  2  3  4  5  6  7  8  9 10 11 12
Original RTs:             1  2  3  4  5  6  7  8  9 10 11 12
Corrected RTs:            1  2  3  4  8  7  6  5  9 10 11 12
@endcode

After correction and sorting, spectra are ordered by actual spotting sequence.

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
