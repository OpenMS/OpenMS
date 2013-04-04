// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>

#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

using namespace OpenMS;
using namespace std;

/**
  @page TOPP_HighResPrecursorMassCorrector HighResPrecursorMassCorrector

  @brief Corrects the precursor mz of high resolution data.

 <CENTER>
 <table>
   <tr>
     <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
     <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ HighResPrecursorMassCorrector \f$ \longrightarrow \f$</td>
     <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
   </tr>
   <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
   </tr>
 </table>
 </CENTER>

  This tool performs precursor mz correction on picked (=centroided) high resolution data.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_HighResPrecursorMassCorrector.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_HighResPrecursorMassCorrector.html
*/

/// @cond TOPPCLASSES

class TOPPHiResPrecursorMassCorrector :
  public TOPPBase
{
public:
  TOPPHiResPrecursorMassCorrector() :
    TOPPBase("HighResPrecursorMassCorrector", "Corrects the precursor mz determined by the instrument software.")
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    // input files
    registerInputFile_("in", "<file>", "", "input file (centroided data)");
    setValidFormats_("in", StringList::create("mzML"));
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", StringList::create("mzML"));
    registerOutputFile_("out_csv", "<file>", "", "Optional csv output file containing columns: precursor rt, uncorrected mz, corrected mz, delta mz\n", false);
    setValidFormats_("out_csv", StringList::create("csv"));
  }

  void getPrecursors_(const PeakMap & exp, vector<Precursor> & precursors, vector<double> & precursors_rt)
  {
    for (Size i = 0; i != exp.size(); ++i)
    {
      vector<Precursor> pcs = exp[i].getPrecursors();
      if (pcs.empty())
      {
        continue;
      }
      vector<double> pcs_rt(pcs.size(), exp[i].getRT());
      copy(pcs.begin(), pcs.end(), back_inserter(precursors));
      copy(pcs_rt.begin(), pcs_rt.end(), back_inserter(precursors_rt));
    }
  }

  void writeHist(String out_csv, const vector<DoubleReal> & deltaMZs, const vector<DoubleReal> & mzs, const vector<DoubleReal> & rts)
  {
    //cout << "writting data" << endl;
    ofstream csv_file(out_csv.c_str());
    csv_file << setprecision(9);

    // header
    csv_file << "RT\tuncorrectedMZ\tcorrectedMZ\tdeltaMZ" << endl;

    // entries
    for (vector<DoubleReal>::const_iterator it = deltaMZs.begin(); it != deltaMZs.end(); ++it)
    {
      UInt index = it - deltaMZs.begin();
      csv_file << rts[index] << "\t" << mzs[index] << "\t" << mzs[index] + *it  << "\t" << *it << endl;
    }
    csv_file.close();
  }

protected:
  void correct(PeakMap & exp, vector<DoubleReal> & deltaMZs, vector<DoubleReal> & mzs, vector<DoubleReal> & rts)
  {
    // load experiment and extract precursors
    vector<Precursor> precursors;  // precursor
    vector<double> precursors_rt;  // RT of precursor MS2 spectrum
    getPrecursors_(exp, precursors, precursors_rt);

    for (Size i = 0; i != precursors_rt.size(); ++i)
    {
      // get precursor rt
      DoubleReal rt = precursors_rt[i];

      // get precursor MZ
      DoubleReal mz = precursors[i].getMZ();

      //cout << rt << " " << mz << endl;

      // get precursor spectrum
      MSExperiment<Peak1D>::ConstIterator rt_it = exp.RTBegin(rt);

      // store index of MS2 spectrum
      UInt precursor_spectrum_idx = rt_it - exp.begin();

      // get parent (MS1) of precursor spectrum
      rt_it = exp.getPrecursorSpectrum(rt_it);

      if (rt_it->getMSLevel() != 1)
      {
        cout << "Error: no MS1 spectrum for this precursor" << endl;
      }

      //cout << rt_it->getRT() << " " << rt_it->size() << endl;

      // find peak (index) closest to expected position
      Size nearest_peak_idx = rt_it->findNearest(mz);

      // get actual position of closest peak
      DoubleReal nearest_peak_mz = (*rt_it)[nearest_peak_idx].getMZ();

      // calculate error between expected and actual position
      DoubleReal nearestPeakError = abs(nearest_peak_mz - mz);

      // check if error is small enough
      if (nearestPeakError < 0.1)
      {
        // sanity check: do we really have the same precursor in the original and the picked spectrum
        if (fabs(exp[precursor_spectrum_idx].getPrecursors()[0].getMZ() - mz) > 0.0001)
        {
          cout << "Error: index is referencing different precursors in original and picked spectrum." << endl;
        }

        // cout << mz << " -> " << nearest_peak_mz << endl;
        DoubleReal deltaMZ = nearest_peak_mz - mz;
        deltaMZs.push_back(deltaMZ);
        mzs.push_back(mz);
        rts.push_back(rt);
        // correct entries
        Precursor corrected_prec = precursors[i];
        corrected_prec.setMZ(nearest_peak_mz);
        exp[precursor_spectrum_idx].getPrecursors()[0] = corrected_prec;
      }
    }
  }

  ExitCodes main_(int, const char **)
  {
    const string in_mzml(getStringOption_("in"));
    const string out_mzml(getStringOption_("out"));
    const string out_csv = getStringOption_("out_csv");

    PeakMap exp;
    MzMLFile().load(in_mzml, exp);

    cout << setprecision(12);

    // determine accuracy
    vector<DoubleReal> deltaMZs;
    vector<DoubleReal> mzs;
    vector<DoubleReal> rts;

    correct(exp, deltaMZs, mzs, rts);

    MzMLFile().store(out_mzml, exp);

    if (out_csv != "")
    {
      writeHist(out_csv, deltaMZs, mzs, rts);
    }

    return EXECUTION_OK;
  }

};

/// @endcond

int main(int argc, const char ** argv)
{
  TOPPHiResPrecursorMassCorrector tool;
  return tool.main(argc, argv);
}
