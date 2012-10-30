// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace OpenMS;

class TOPPRNPxlXICFilter :
  public TOPPBase
{
public:
  TOPPRNPxlXICFilter() :
    TOPPBase("RNPxlXICFilter", "Remove MS2 spectra from treatment based on the fold change between control and treatment.", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    // input files
    registerInputFile_("control", "<file>", "", "input mzML file");
    registerInputFile_("treatment", "<file>", "", "input mzML file");
    registerDoubleOption_("fold_change", "", 2.0, "fold change between XICs", false, false);
    registerDoubleOption_("rt_tol", "", 20, "RT tolerance in [s] for finding max peak (whole RT range around RT middle)", false, false);
    registerDoubleOption_("mz_tol", "", 10, "m/z tolerance in [ppm] for finding a peak", false, false);

    // output files
    registerOutputFile_("out", "<file>", "", "output file");
  }

  void filterByFoldChange(const MSExperiment<>& exp1, const MSExperiment<>& exp2,
                          const vector<DoubleReal>& pc_ms2_rts, const vector<DoubleReal>& pc_mzs,
                          const DoubleReal rttol, const DoubleReal mztol, DoubleReal fold_change,
                          vector<DoubleReal>& control_XIC_larger,
                          vector<DoubleReal>& treatment_XIC_larger,
                          vector<DoubleReal>& indifferent_XICs)
  {
    assert(pc_mzs.size() == pc_ms2_rts.size());

    // search for each EIC and add up
    for (Size i = 0; i < pc_mzs.size(); ++i)
    {
      //cerr << "start" << endl;
      DoubleReal pc_ms2_rt = pc_ms2_rts[i];
      DoubleReal pc_mz = pc_mzs[i];

      //std::cerr << "Rt" << cm[i].getRT() << "  mz: " << cm[i].getMZ() << " R " <<  cm[i].getMetaValue("rank") << "\n";

      DoubleReal mz_da = mztol * pc_mzs[i] / 1e6; // mz tolerance in Dalton
      DoubleReal rt_start = pc_ms2_rts[i] - rttol / 2.0;

      // get area iterator (is MS1 only!) for rt and mz window
      MSExperiment<>::ConstAreaIterator it1 = exp1.areaBeginConst(pc_ms2_rt - rttol / 2, pc_ms2_rt + rttol / 2, pc_mz - mz_da, pc_mz  + mz_da);
      MSExperiment<>::ConstAreaIterator it2 = exp2.areaBeginConst(pc_ms2_rt - rttol / 2, pc_ms2_rt + rttol / 2, pc_mz - mz_da, pc_mz  + mz_da);

      // determine maximum number of MS1 scans in retention time window
      set<DoubleReal> rts1;
      set<DoubleReal> rts2;
      for (; it1 != exp1.areaEndConst(); ++it1)
      {
        rts1.insert(it1.getRT());
      }

      for (; it2 != exp2.areaEndConst(); ++it2)
      {
        rts2.insert(it2.getRT());
      }

      Size length = std::max(rts1.size(), rts2.size()) / 2.0;

      cout << length << endl;
      if (length == 0)
      {
        cerr << "WARNING: no MS1 scans in retention time window found in both maps (mz: " << pc_mzs[i] << " / rt: " << pc_ms2_rts[i] << ")" << endl;
        continue;
      }

      vector<DoubleReal> XIC1(length, 0.0);
      vector<DoubleReal> XIC2(length, 0.0);

      it1 = exp1.areaBeginConst(pc_ms2_rt - rttol / 2, pc_ms2_rt + rttol / 2, pc_mz - mz_da, pc_mz + mz_da);
      it2 = exp2.areaBeginConst(pc_ms2_rt - rttol / 2, pc_ms2_rt + rttol / 2, pc_mz - mz_da, pc_mz + mz_da);

      for (; it1 != exp1.areaEndConst(); ++it1)
      {
        DoubleReal relative_rt = (it1.getRT() - rt_start) / rttol;
        Size bin = relative_rt * (length - 1);
        XIC1[bin] += it1->getIntensity();
        if (bin >= length)
        {
          bin = length - 1;
        }

      }

      for (; it2 != exp2.areaEndConst(); ++it2)
      {
        DoubleReal relative_rt = (it2.getRT() - rt_start) / rttol;
        Size bin = relative_rt * (length - 1);
        if (bin >= length)
        {
          bin = length - 1;
        }
        XIC2[bin] += it2->getIntensity();
      }

      DoubleReal total_itensity1 = std::accumulate(XIC1.begin(), XIC1.end(), 0.0);
      DoubleReal total_itensity2 = std::accumulate(XIC2.begin(), XIC2.end(), 0.0);

      DoubleReal ratio = total_itensity2 / (total_itensity1 + 1);

      //cout << pc_ms2_rt << "/" << pc_mz << " has ratio: " << ratio << " determined on " << length << " bins" << endl;

      if (ratio < 1.0 / fold_change)
      {
        control_XIC_larger.push_back(pc_ms2_rt);
      }
      else if (ratio > fold_change)
      {
        treatment_XIC_larger.push_back(pc_ms2_rt);
      }
      else
      {
        indifferent_XICs.push_back(pc_ms2_rt);
        continue;
      }

      for (Size k = 0; k != length; ++k)
      {
        cout << k << ": " << rt_start + rttol / length * k  << ": " << XIC1[k] << " " << XIC2[k] << endl;
      }
    }

    cout << "control larger: " << control_XIC_larger.size() << " treatment larger: " << treatment_XIC_larger.size() << " indifferent: " << indifferent_XICs.size() << endl;

    return;
  }

  ExitCodes main_(int, const char**)
  {
    // Parameter parsing
    const string control_mzml(getStringOption_("control"));
    const string treatment_mzml(getStringOption_("treatment"));
    const string out_mzml(getStringOption_("out"));
    const DoubleReal mz_tolerance_ppm = getDoubleOption_("mz_tol");
    const DoubleReal fold_change = getDoubleOption_("fold_change");
    const DoubleReal rt_tolerance_s = getDoubleOption_("rt_tol");

    // load experiments
    MSExperiment<> exp_control;
    MzMLFile mzml_file;
    mzml_file.load(control_mzml, exp_control);

    MSExperiment<> exp_treatment;
    mzml_file.load(treatment_mzml, exp_treatment);

    // extract precursor mz and rts
    vector<DoubleReal> pc_mzs;
    vector<DoubleReal> pc_ms2_rts;
    for (Size i = 0; i != exp_treatment.size(); ++i)
    {
      if (exp_treatment[i].getMSLevel() == 2)
      {
        if (!exp_treatment[i].getPrecursors().empty())
        {
          // cout << i << endl;
          DoubleReal pc_mz = exp_treatment[i].getPrecursors()[0].getMZ();
          DoubleReal ms2_rt = exp_treatment[i].getRT(); // use rt of MS2
          pc_mzs.push_back(pc_mz);
          pc_ms2_rts.push_back(ms2_rt);
        }
      }
    }

    vector<DoubleReal> control_XIC_larger_rts;
    vector<DoubleReal> treatment_XIC_larger_rts;
    vector<DoubleReal> indifferent_XICs_rts;

    filterByFoldChange(exp_control, exp_treatment,
                       pc_ms2_rts, pc_mzs,
                       rt_tolerance_s, mz_tolerance_ppm, fold_change,
                       control_XIC_larger_rts,
                       treatment_XIC_larger_rts,
                       indifferent_XICs_rts);


    MSExperiment<> exp_out = exp_treatment;
    exp_out.clear(false); // don't clear meta-data

    for (Size i = 0; i != exp_treatment.size(); ++i)
    {
      Size ms_level = exp_treatment[i].getMSLevel();

      if (ms_level == 1)
      {
        exp_out.push_back(exp_treatment[i]);
        continue;
      }
      else if (ms_level == 2)
      {
        // determine if pc is in list -> passed
        DoubleReal rt = exp_treatment[i].getRT();
        for (Size j = 0; j != treatment_XIC_larger_rts.size(); ++j)
        {
          if (fabs(rt - treatment_XIC_larger_rts[j]) <= 0.001)
          {
            DoubleReal pc_mz = exp_treatment[i].getPrecursors()[0].getMZ();
            DoubleReal pc_charge = exp_treatment[i].getPrecursors()[0].getCharge();
            DoubleReal pc_mass = pc_mz * pc_charge - pc_charge * Constants::PROTON_MASS_U;

            exp_out.push_back(exp_treatment[i]);
            break;
          }
        }
      }
    }

    mzml_file.store(out_mzml, exp_out);

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPRNPxlXICFilter tool;
  return tool.main(argc, argv);
}
