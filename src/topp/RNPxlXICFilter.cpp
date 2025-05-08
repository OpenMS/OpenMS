// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
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

/**
@page TOPP_RNPxlXICFilter RNPxlXICFilter
@brief Filters MS2 spectra based on XIC intensities in control and treatment. Used in RNPxl experiments to reduce candidate spectra.

<CENTER>
<table>
<tr>
  <th ALIGN = "center"> potential predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=2> &rarr; RNPxlXICFilter &rarr;</td>
  <th ALIGN = "center"> potential successor tools </td>
</tr>
<tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> NuXL </td>
</tr>
</table>
</CENTER>

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_RNPxlXICFilter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_RNPxlXICFilter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPRNPxlXICFilter :
  public TOPPBase
{
public:
  TOPPRNPxlXICFilter() :
    TOPPBase("RNPxlXICFilter", "Remove MS2 spectra from treatment based on the fold change between control and treatment.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    // input files

    registerInputFile_("control", "<file>", "", "input mzML file");
    setValidFormats_("control", ListUtils::create<String>("mzML"));
    registerInputFile_("treatment", "<file>", "", "input mzML file");
    setValidFormats_("treatment", ListUtils::create<String>("mzML"));

    registerDoubleOption_("fold_change", "", 2.0, "fold change between XICs", false, false);
    registerDoubleOption_("rt_tol", "", 20, "RT tolerance in [s] for finding max peak (whole RT range around RT middle)", false, false);
    registerDoubleOption_("mz_tol", "", 10, "m/z tolerance in [ppm] for finding a peak", false, false);

    // output files
    registerOutputFile_("out", "<file>", "", "output of the treatment file after XIC filtering.");
    setValidFormats_("out", ListUtils::create<String>("mzML"));    
  }

  void filterByFoldChange(const PeakMap& exp1, const PeakMap& exp2,
                          const vector<double>& pc_ms2_rts, const vector<double>& pc_mzs,
                          const double rttol, const double mztol, double fold_change,
                          vector<double>& control_XIC_larger,
                          vector<double>& treatment_XIC_larger,
                          vector<double>& indifferent_XICs)
  {
    assert(pc_mzs.size() == pc_ms2_rts.size());

    // search for each EIC and add up
    for (Size i = 0; i < pc_mzs.size(); ++i)
    {
      //cerr << "start" << endl;
      double pc_ms2_rt = pc_ms2_rts[i];
      double pc_mz = pc_mzs[i];

      //std::cerr << "Rt" << cm[i].getRT() << "  mz: " << cm[i].getMZ() << " R " <<  cm[i].getMetaValue("rank") << "\n";

      double mz_da = mztol * pc_mzs[i] / 1e6; // mz tolerance in Dalton
      double rt_start = pc_ms2_rts[i] - rttol / 2.0;

      // get area iterator (is MS1 only!) for rt and mz window
      PeakMap::ConstAreaIterator it1 = exp1.areaBeginConst(pc_ms2_rt - rttol / 2, pc_ms2_rt + rttol / 2, pc_mz - mz_da, pc_mz  + mz_da);
      PeakMap::ConstAreaIterator it2 = exp2.areaBeginConst(pc_ms2_rt - rttol / 2, pc_ms2_rt + rttol / 2, pc_mz - mz_da, pc_mz  + mz_da);

      // determine maximum number of MS1 scans in retention time window
      set<double> rts1;
      set<double> rts2;
      for (; it1 != exp1.areaEndConst(); ++it1)
      {
        rts1.insert(it1.getRT());
      }

      for (; it2 != exp2.areaEndConst(); ++it2)
      {
        rts2.insert(it2.getRT());
      }

      Size length = std::max(rts1.size(), rts2.size()) / 2.0;

      //cout << length << endl;
      if (length == 0)
      {
        cerr << "WARNING: no MS1 scans in retention time window found in both maps (mz: " << pc_mzs[i] << " / rt: " << pc_ms2_rts[i] << ")" << endl;
        continue;
      }

      vector<double> XIC1(length, 0.0);
      vector<double> XIC2(length, 0.0);

      it1 = exp1.areaBeginConst(pc_ms2_rt - rttol / 2, pc_ms2_rt + rttol / 2, pc_mz - mz_da, pc_mz + mz_da);
      it2 = exp2.areaBeginConst(pc_ms2_rt - rttol / 2, pc_ms2_rt + rttol / 2, pc_mz - mz_da, pc_mz + mz_da);

      for (; it1 != exp1.areaEndConst(); ++it1)
      {
        double relative_rt = (it1.getRT() - rt_start) / rttol;
        Size bin = relative_rt * (length - 1);
        XIC1[bin] += it1->getIntensity();
        if (bin >= length)
        {
          bin = length - 1;
        }

      }

      for (; it2 != exp2.areaEndConst(); ++it2)
      {
        double relative_rt = (it2.getRT() - rt_start) / rttol;
        Size bin = relative_rt * (length - 1);
        if (bin >= length)
        {
          bin = length - 1;
        }
        XIC2[bin] += it2->getIntensity();
      }

      double total_itensity1 = std::accumulate(XIC1.begin(), XIC1.end(), 0.0);
      double total_itensity2 = std::accumulate(XIC2.begin(), XIC2.end(), 0.0);

      double ratio = total_itensity2 / (total_itensity1 + 1);

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
      /*
      for (Size k = 0; k != length; ++k)
      {
        cout << k << ": " << rt_start + rttol / length * k  << ": " << XIC1[k] << " " << XIC2[k] << endl;
      }
      */
    }

    cout << "control larger: " << control_XIC_larger.size() << " treatment larger: " << treatment_XIC_larger.size() << " indifferent: " << indifferent_XICs.size() << endl;

    return;
  }

  ExitCodes main_(int, const char**) override
  {
    // Parameter parsing
    const string control_mzml(getStringOption_("control"));
    const string treatment_mzml(getStringOption_("treatment"));

    const string out_mzml(getStringOption_("out"));
    const double mz_tolerance_ppm = getDoubleOption_("mz_tol");
    const double fold_change = getDoubleOption_("fold_change");
    const double rt_tolerance_s = getDoubleOption_("rt_tol");

    // load experiments
    PeakMap exp_control;
    FileHandler().loadExperiment(control_mzml, exp_control, {FileTypes::MZML});

    PeakMap exp_treatment;
    FileHandler().loadExperiment(treatment_mzml, exp_treatment, {FileTypes::MZML});

    // extract precursor mz and rts
    vector<double> pc_mzs;
    vector<double> pc_ms2_rts;
    for (Size i = 0; i != exp_treatment.size(); ++i)
    {
      if (exp_treatment[i].getMSLevel() == 2)
      {
        if (!exp_treatment[i].getPrecursors().empty())
        {
          // cout << i << endl;
          double pc_mz = exp_treatment[i].getPrecursors()[0].getMZ();
          double ms2_rt = exp_treatment[i].getRT(); // use rt of MS2
          pc_mzs.push_back(pc_mz);
          pc_ms2_rts.push_back(ms2_rt);
        }
      }
    }

    vector<double> control_XIC_larger_rts;
    vector<double> treatment_XIC_larger_rts;
    vector<double> indifferent_XICs_rts;

    filterByFoldChange(exp_control, exp_treatment,
                       pc_ms2_rts, pc_mzs,
                       rt_tolerance_s, mz_tolerance_ppm, fold_change,
                       control_XIC_larger_rts,
                       treatment_XIC_larger_rts,
                       indifferent_XICs_rts);


    PeakMap exp_out = exp_treatment;
    exp_out.clear(false); // don't clear meta-data

    for (Size i = 0; i != exp_treatment.size(); ++i)
    {
      Size ms_level = exp_treatment[i].getMSLevel();

      if (ms_level == 1)
      {
        exp_out.addSpectrum(exp_treatment[i]);
        continue;
      }
      else if (ms_level == 2)
      {
        // determine if pc is in list -> passed
        double rt = exp_treatment[i].getRT();
        for (Size j = 0; j != treatment_XIC_larger_rts.size(); ++j)
        {
          if (fabs(rt - treatment_XIC_larger_rts[j]) <= 0.001)
          {
            exp_out.addSpectrum(exp_treatment[i]);
            break;
          }
        }
      }
    }

    FileHandler().storeExperiment(out_mzml, exp_out, {FileTypes::MZML});

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPRNPxlXICFilter tool;
  return tool.main(argc, argv);
}

///@endcond
