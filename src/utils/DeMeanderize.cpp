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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_DeMeanderize DeMeanderize

    @brief Repairs MALDI experiments which were spotted line by line.



    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_DeMeanderize.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_DeMeanderize.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDeMeanderize :
  public TOPPBase
{
public:
  TOPPDeMeanderize() :
    TOPPBase("DeMeanderize", "Orders the spectra of MALDI spotting plates correctly.", false)
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
    MzMLFile f;
    f.setLogType(log_type_);
    f.load(in, exp);

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


    f.store(out, exp);

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPDeMeanderize tool;
  return tool.main(argc, argv);
}

/// @endcond
