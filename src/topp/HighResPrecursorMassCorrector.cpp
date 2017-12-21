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
// $Authors: Timo Sachsenberg, Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
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
     <td VALIGN= "middle" ROWSPAN=2> \f$ \longrightarrow \f$ HighResPrecursorMassCorrector \f$ \longrightarrow \f$</td>
     <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
   </tr>
   <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
   </tr>
 </table>
 </CENTER>

  This tool performs precursor m/z correction on picked (=centroided) high resolution data.

  Three methods are available: 'nearest_peak', 'highest_intensity_peak' and 'feature'.
    - nearest_peak: Use nearest centroided MS1 peak for precursor mass correction.
    - highest_intensity_peak: Use highest intensity centroided MS1 peak in a given mass range for precursor mass correction.
    - feature: Use features for precursor mass correction, which allows for charge correction.

  The method hightest_intensity_peak searches in a specific m/z-window of the precursor information for the peak with the highest intensity.
  Suggestioned value 1/maximal expected charge. E.g maximal expected charge 5, m/z-window = +/- 0.2 Da

  See the corresponding parameter subsection for details.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_HighResPrecursorMassCorrector.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_HighResPrecursorMassCorrector.html
*/

/// @cond TOPPCLASSES

const string csv_header = "RT,uncorrectedMZ,correctedMZ,deltaMZ";

class TOPPHiResPrecursorMassCorrector :
    public TOPPBase
{
  public:
    TOPPHiResPrecursorMassCorrector() :
      TOPPBase("HighResPrecursorMassCorrector", "Corrects the precursor mass and charge determined by the instrument software.")
    {
    }

  protected:
    void registerOptionsAndFlags_() override
    {
      // input files
      registerInputFile_("in", "<file>", "", "Input file (centroided data)");
      setValidFormats_("in", ListUtils::create<String>("mzML"));

      registerOutputFile_("out", "<file>", "", "Output file");
      setValidFormats_("out", ListUtils::create<String>("mzML"));
      
      registerTOPPSubsection_("feature", "Use features for precursor mass correction.");
      registerInputFile_("feature:in", "<file>", "", "Features used to correct precursor masses.", false);
      setValidFormats_("feature:in", ListUtils::create<String>("featureXML"));
      registerDoubleOption_("feature:mz_tolerance", "<num>", 5.0, "The precursor mass tolerance. Used to determine matching to feature mass traces.", false);
      registerStringOption_("feature:mz_tolerance_unit", "<choice>", "ppm", "Unit of precursor mass tolerance", false);
      setValidStrings_("feature:mz_tolerance_unit", ListUtils::create<String>("Da,ppm"));
      registerDoubleOption_("feature:rt_tolerance", "<num>", 0.0, "Additional retention time tolerance added to feature boundaries.", false);
      registerIntOption_("feature:max_trace", "<num>", 2, "Maximum isotopic trace considered in matching a precursor to a feature.", false, true);
      registerFlag_("feature:believe_charge", "Assume precursor charge to be correct.");
      registerFlag_("feature:keep_original", "Make a copy of the precursor and MS2 (true) or discard the original (false).");
      registerFlag_("feature:assign_all_matching", "Correct a precursor using all matching features (true) or only the nearest (false). Only evaluated if copies are created (feature:keep_original).");

      registerTOPPSubsection_("nearest_peak", "Use nearest centroided MS1 peak for precursor mass correction.");
      registerDoubleOption_("nearest_peak:mz_tolerance", "<num>", 0.0, "The precursor mass tolerance to find the closest MS1 peak. (Disable method by setting value to 0.0)", false);
      registerStringOption_("nearest_peak:mz_tolerance_unit", "<choice>", "ppm", "Unit of precursor mass tolerance", false);
      setValidStrings_("nearest_peak:mz_tolerance_unit", ListUtils::create<String>("Da,ppm"));

      registerTOPPSubsection_("highest_intensity_peak", "Use centroided MS1 peak with the highest intensity in a certrain mass range - for precursor mass correction");
      registerDoubleOption_("highest_intensity_peak:mz_tolerance", "<num>", 0.0, "The precursor mass tolerance to find the highest intensity MS1 peak (Da). Suggested value 1/max. expected charge. (Disable method by setting value to 0.0)", false);

      registerOutputFile_("out_csv", "<file>", "", "Optional CSV output file for results on 'nearest_peak' or 'highest_intensity_peak' algorithm (see corresponding subsection) containing columns: " + ListUtils::concatenate(ListUtils::create<String>(csv_header), ", ") + ".", false);
      setValidFormats_("out_csv", ListUtils::create<String>("csv"));
    }

    void getPrecursors_(const PeakMap & exp, vector<Precursor> & precursors, vector<double> & precursors_rt, vector<Size> precursor_scan_index)
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
        precursor_scan_index.push_back(i);
      }
    }

    void writeHist(const String& out_csv, const vector<double> & deltaMZs, const vector<double> & mzs, const vector<double> & rts)
    {
      //cout << "writting data" << endl;
      ofstream csv_file(out_csv.c_str());
      csv_file << setprecision(9);

      // header
      csv_file << ListUtils::concatenate(ListUtils::create<String>(csv_header), "\t") << "\n";

      // entries
      for (vector<double>::const_iterator it = deltaMZs.begin(); it != deltaMZs.end(); ++it)
      {
        UInt index = it - deltaMZs.begin();
        csv_file << rts[index] << "\t" << mzs[index] << "\t" << mzs[index] + *it  << "\t" << *it << "\n";
      }
      csv_file.close();
    }

  protected:
    bool overlaps_(const Feature& feature, const double rt, const double pc_mz, const double rt_tolerance) const
    {
      if (feature.getConvexHulls().empty())
      {
        LOG_WARN << "HighResPrecursorMassCorrector warning: at least one feature has no convex hull - omitting feature for matching" << std::endl;
      }

      // get bounding box and extend by retention time tolerance
      DBoundingBox<2> box = feature.getConvexHull().getBoundingBox();
      DPosition<2> extend_rt(rt_tolerance, 0.01);
      box.setMin(box.minPosition() - extend_rt);
      box.setMax(box.maxPosition() + extend_rt);

      DPosition<2> pc_pos(rt, pc_mz);
      if (box.encloses(pc_pos))
      {
        return true;
      }
      else
      {
        return false;
      }
    }

    bool compatible_(const Feature& feature, double pc_mz, double mz_tolerance, Size max_trace_number = 2)
    {
      const int f_charge = feature.getCharge();
      const double f_mz = feature.getMZ();
      double trace = Math::round((pc_mz - f_mz) / (Constants::C13C12_MASSDIFF_U / f_charge)); // isotopic trace number at precursor mz
      double mass_error = fabs(pc_mz - (f_mz + trace * (Constants::C13C12_MASSDIFF_U / f_charge)));

      if (mass_error < mz_tolerance && (trace < max_trace_number + 0.01))
      {
        if (debug_level_ > 1)
        {
          LOG_INFO << "trace: " << (int)(trace + 0.5) << " feature_rt:" << feature.getRT() << " feature_mz:" << feature.getMZ() << " precursor_mz:" << pc_mz << endl;
        }
        return true;
      }
      else
      {
        return false;
      }
    }

    set<Size> correctToNearestMS1Peak(PeakMap & exp, double mz_tolerance, bool ppm, vector<double> & deltaMZs, vector<double> & mzs, vector<double> & rts)
    {
      set<Size> corrected_precursors;
      // load experiment and extract precursors
      vector<Precursor> precursors;  // precursor
      vector<double> precursors_rt;  // RT of precursor MS2 spectrum
      vector<Size> precursor_scan_index;
      getPrecursors_(exp, precursors, precursors_rt, precursor_scan_index);

      for (Size i = 0; i != precursors_rt.size(); ++i)
      {
        // get precursor rt
        double rt = precursors_rt[i];

        // get precursor MZ
        double mz = precursors[i].getMZ();

        //cout << rt << " " << mz << endl;

        // get precursor spectrum
        PeakMap::ConstIterator rt_it = exp.RTBegin(rt - 1e-8);

        // store index of MS2 spectrum
        UInt precursor_spectrum_idx = rt_it - exp.begin();

        // get parent (MS1) of precursor spectrum
        rt_it = exp.getPrecursorSpectrum(rt_it);

        if (rt_it->getMSLevel() != 1)
        {
          LOG_WARN << "Error: no MS1 spectrum for this precursor" << endl;
        }

        //cout << rt_it->getRT() << " " << rt_it->size() << endl;

        // find peak (index) closest to expected position
        Size nearest_peak_idx = rt_it->findNearest(mz);

        // get actual position of closest peak
        double nearest_peak_mz = (*rt_it)[nearest_peak_idx].getMZ();

        // calculate error between expected and actual position
        double nearestPeakError = ppm ? abs(nearest_peak_mz - mz)/mz * 1e6 : abs(nearest_peak_mz - mz);

        // check if error is small enough
        if (nearestPeakError < mz_tolerance)
        {
          // sanity check: do we really have the same precursor in the original and the picked spectrum
          if (fabs(exp[precursor_spectrum_idx].getPrecursors()[0].getMZ() - mz) > 0.0001)
          {
            LOG_WARN << "Error: index is referencing different precursors in original and picked spectrum." << endl;
          }

          // cout << mz << " -> " << nearest_peak_mz << endl;
          double deltaMZ = nearest_peak_mz - mz;
          deltaMZs.push_back(deltaMZ);
          mzs.push_back(mz);
          rts.push_back(rt);
          // correct entries
          Precursor corrected_prec = precursors[i];
          corrected_prec.setMZ(nearest_peak_mz);
          exp[precursor_spectrum_idx].getPrecursors()[0] = corrected_prec;
          corrected_precursors.insert(precursor_spectrum_idx);
        }
      }
      return corrected_precursors;
    }

    //Selection of the peak with the highest intensity as corrected precursor mass in a given mass range (e.g. precursor mass +/- 0.2 Da)
    set<Size> correctToHighestIntensityMS1Peak(PeakMap & exp, double mz_tolerance, vector<double> & deltaMZs, vector<double> & mzs, vector<double> & rts)
    {
      set<Size> corrected_precursors;
      // load experiment and extract precursors
      vector<Precursor> precursors;  // precursor
      vector<double> precursors_rt;  // RT of precursor MS2 spectrum
      vector<Size> precursor_scan_index;
      getPrecursors_(exp, precursors, precursors_rt, precursor_scan_index);
      int count_error_highest_intenstiy = 0;

      for (Size i = 0; i != precursors_rt.size(); ++i)
      {
        // get precursor rt
        double rt = precursors_rt[i];

        // get precursor MZ
        double mz = precursors[i].getMZ();

        // cout << rt << " " << mz << endl;

        // retrieves iterator of the MS2 fragment sprectrum
        PeakMap::ConstIterator rt_it = exp.RTBegin(rt - 1e-8);

        // store index of MS2 spectrum
        UInt precursor_spectrum_idx = rt_it - exp.begin();

        // get parent (MS1) of precursor spectrum
        rt_it = exp.getPrecursorSpectrum(rt_it);

        if (rt_it->getMSLevel() != 1)
        {
          LOG_WARN << "Error: no MS1 spectrum for this precursor" << endl;
        }

        MSSpectrum::ConstIterator left = rt_it->MZBegin(mz - mz_tolerance);
        MSSpectrum::ConstIterator right = rt_it->MZEnd(mz + mz_tolerance);

        // no MS1 precursor peak in +- tolerance window found
        if (left == right || left->getMZ() > mz + mz_tolerance)
        {
          count_error_highest_intenstiy += 1;
        }

        MSSpectrum::ConstIterator max_intensity_it = std::max_element(left, right, Peak1D::IntensityLess());

        // find peak (index) with highest intensity to expected position
        Size highest_peak_idx = max_intensity_it - rt_it->begin();

        // get actual position of highest intensity peak
        double highest_peak_mz = (*rt_it)[highest_peak_idx].getMZ();

        // cout << mz << " -> " << nearest_peak_mz << endl;
        double deltaMZ = highest_peak_mz - mz;
        deltaMZs.push_back(deltaMZ);
        mzs.push_back(mz);
        rts.push_back(rt);
        // correct entries
        Precursor corrected_prec = precursors[i];
        corrected_prec.setMZ(highest_peak_mz);
        exp[precursor_spectrum_idx].getPrecursors()[0] = corrected_prec;
        corrected_precursors.insert(precursor_spectrum_idx);

      }

      if (count_error_highest_intenstiy != 0)
      {
        LOG_WARN << "Error: The method highest_intensity_peak failed" << count_error_highest_intenstiy << "times.";
      }

      return corrected_precursors;
    }

    // Wrong assignment of the mono-isotopic mass for precursors are assumed:
    // - if precursor_mz matches the mz of a non-monoisotopic feature mass trace
    // - and in the case that believe_charge is true: if feature_charge matches the precursor_charge
    // In the case of wrong mono-isotopic assignment several options for correction are available:
    // keep_original will create a copy of the precursor and tandem spectrum for the new mono-isotopic mass trace and retain the original one
    // all_matching_features does this not for only the closest feature but all features in a question
    set<Size> correctToNearestFeature(const FeatureMap& features, PeakMap & exp, double rt_tolerance_s = 0.0, double mz_tolerance = 0.0, bool ppm = true, bool believe_charge = false, bool keep_original = false, bool all_matching_features = false, int max_trace = 2)
    {
      set<Size> corrected_precursors;
      // for each precursor/MS2 find all features that are in the given tolerance window (bounding box + rt tolerances)
      // if believe_charge is set, only add features that match the precursor charge
      map<Size, set<Size> > scan_idx_to_feature_idx;

      for (Size scan = 0; scan != exp.size(); ++scan)
      {
        // skip non-tandem mass spectra
        if (exp[scan].getMSLevel() != 2 || exp[scan].getPrecursors().empty()) continue;

        // extract precusor / MS2 information
        const double pc_mz = exp[scan].getPrecursors()[0].getMZ();
        const double rt = exp[scan].getRT();
        const int pc_charge = exp[scan].getPrecursors()[0].getCharge();

        for (Size f = 0; f != features.size(); ++f)
        {
          // feature  is incompatible if believe_charge is set and charges don't match
          if (believe_charge && features[f].getCharge() != pc_charge) continue;

          // check if precursor/MS2 position overlap with feature
          if (overlaps_(features[f], rt, pc_mz, rt_tolerance_s))
          {
            scan_idx_to_feature_idx[scan].insert(f);
          }
        }
      }

      // filter sets to retain compatible features:
      // if precursor_mz = feature_mz + n * feature_charge (+/- mz_tolerance) a feature is compatible, others are removed from the set
      for (map<Size, set<Size> >::iterator it = scan_idx_to_feature_idx.begin(); it != scan_idx_to_feature_idx.end(); ++it)
      {
        const Size scan = it->first;
        const double pc_mz = exp[scan].getPrecursors()[0].getMZ();
        const double mz_tolerance_da = ppm ? pc_mz * mz_tolerance * 1e-6  : mz_tolerance;

        // Note: This is the "delete while iterating" pattern so mind the pre- and postincrement
        for (set<Size>::iterator sit = it->second.begin(); sit != it->second.end(); )
        {
          if (!compatible_(features[*sit], pc_mz, mz_tolerance_da, max_trace))
          {
            it->second.erase(sit++);
          }
          else
          {
            ++sit;
          }
        }
      }

      // remove entries with no compatible features (empty sets).
      // Note: This is the "delete while iterating" pattern so mind the pre- and postincrement
      for (map<Size, set<Size> >::iterator it = scan_idx_to_feature_idx.begin(); it != scan_idx_to_feature_idx.end(); )
      {
        if (it->second.empty())
        {
          scan_idx_to_feature_idx.erase(it++);
        }
        else
        {
          ++it;
        }
      }

      if (debug_level_ > 0)
      {
        LOG_INFO << "Number of precursors with compatible features: " << scan_idx_to_feature_idx.size() << endl;
      }

      if (!all_matching_features)
      {
        // keep only nearest features in set
        for (map<Size, set<Size> >::iterator it = scan_idx_to_feature_idx.begin(); it != scan_idx_to_feature_idx.end(); ++it)
        {
          const Size scan = it->first;
          const double pc_rt = exp[scan].getRT();

          double min_distance = 1e16;
          set<Size>::iterator best_feature = it->second.begin();

          // determine nearest/best feature
          for (set<Size>::iterator sit = it->second.begin(); sit != it->second.end(); ++sit)
          {
            const double current_distance = fabs(pc_rt - features[*sit].getRT());
            if (current_distance < min_distance)
            {
              min_distance = current_distance;
              best_feature = sit;
            }
          }

          // delete all except the nearest/best feature
          // Note: This is the "delete while iterating" pattern so mind the pre- and postincrement
          for (set<Size>::iterator sit = it->second.begin(); sit != it->second.end(); )
          {
            if (sit != best_feature)
            {
              it->second.erase(sit++);
            }
            else
            {
              ++sit;
            }
          }
        }
      }

      // depending on all_matching_features option, only the nearest or all features are contained in the sets
      // depending on options: move/copy corrected precursor and tandem spectrum
      if (keep_original)
      {
        // duplicate spectra for each feature in set and adapt precursor_mz and precursor_charge to feature_mz and feature_charge
        for (map<Size, set<Size> >::iterator it = scan_idx_to_feature_idx.begin(); it != scan_idx_to_feature_idx.end(); ++it)
        {
          const Size scan = it->first;
          MSSpectrum spectrum = exp[scan];
          corrected_precursors.insert(scan);
          for (set<Size>::iterator f_it = it->second.begin(); f_it != it->second.end(); ++f_it)
          {
            spectrum.getPrecursors()[0].setMZ(features[*f_it].getMZ());
            spectrum.getPrecursors()[0].setCharge(features[*f_it].getCharge());
            exp.addSpectrum(spectrum);
          }
        }
      }
      else
      {
        // set precursor_mz and _charge to the feature_mz and _charge
        for (map<Size, set<Size> >::iterator it = scan_idx_to_feature_idx.begin(); it != scan_idx_to_feature_idx.end(); ++it)
        {
          const Size scan = it->first;
          exp[scan].getPrecursors()[0].setMZ(features[*it->second.begin()].getMZ());
          exp[scan].getPrecursors()[0].setCharge(features[*it->second.begin()].getCharge());
          corrected_precursors.insert(scan);
        }
      }
      return corrected_precursors;
    }

    ExitCodes main_(int, const char **) override
    {
      const string in_mzml(getStringOption_("in"));
      const string in_feature(getStringOption_("feature:in"));
      const string out_mzml(getStringOption_("out"));
      const string out_csv = getStringOption_("out_csv");

      const double mz_tolerance = getDoubleOption_("feature:mz_tolerance");
      const bool mz_unit_ppm = getStringOption_("feature:mz_tolerance_unit") == "ppm" ? true : false;
      const double rt_tolerance = getDoubleOption_("feature:rt_tolerance");
      const int max_trace = getIntOption_("feature:max_trace");
      bool keep_original = getFlag_("feature:keep_original");
      bool assign_all_matching = getFlag_("feature:assign_all_matching");
      bool believe_charge = getFlag_("feature:believe_charge");

      const double nearest_peak_mz_tolerance = getDoubleOption_("nearest_peak:mz_tolerance");
      const bool nearest_peak_ppm = getStringOption_("nearest_peak:mz_tolerance_unit") == "ppm" ? true : false;

      const double highest_intensity_peak_mz_tolerance = getDoubleOption_("highest_intensity_peak:mz_tolerance");

      PeakMap exp;
      MzMLFile().load(in_mzml, exp);

      cout << setprecision(12);

      // determine accuracy
      vector<double> deltaMZs;
      vector<double> mzs;
      vector<double> rts;
      set<Size> corrected_precursors; // spectrum index of corrected precursors

      if ((nearest_peak_mz_tolerance <= 0.0) && (highest_intensity_peak_mz_tolerance <= 0.0) && in_feature.empty())
      {
        LOG_ERROR << "No method for PC correction requested. Either provide featureXML input files or set 'nearest_peak:mz_tolerance' > 0 or specify a 'highest_intensity_peak:mz_tolerance' > 0" << std::endl;
        return MISSING_PARAMETERS;
      }

      // perform correction to closest MS1 peak
      set<Size> corrected_to_nearest_peak;
      if (nearest_peak_mz_tolerance > 0.0 && highest_intensity_peak_mz_tolerance <= 0.0)
      {
        corrected_to_nearest_peak = correctToNearestMS1Peak(exp, nearest_peak_mz_tolerance, nearest_peak_ppm, deltaMZs, mzs, rts);
      }

      //perform correction to highest intensity MS1 peak
      set<Size> corrected_to_highest_intensity_peak;
      if (highest_intensity_peak_mz_tolerance > 0.0)
      {
        corrected_to_highest_intensity_peak = correctToHighestIntensityMS1Peak(exp, highest_intensity_peak_mz_tolerance, deltaMZs, mzs, rts);
      }
 
      // perform correction to closest feature (also corrects charge if not disabled)
      set<Size> corrected_to_nearest_feature;      
      if (!in_feature.empty())
      {
        FeatureMap features;
        FeatureXMLFile().load(in_feature, features);
        corrected_to_nearest_feature = correctToNearestFeature(features, exp, rt_tolerance, mz_tolerance, mz_unit_ppm, believe_charge, keep_original, assign_all_matching, max_trace);
        corrected_precursors.insert(corrected_to_nearest_feature.begin(), corrected_to_nearest_feature.end());
      }

      MzMLFile().store(out_mzml, exp);

      if (!out_csv.empty())
      {
        if (nearest_peak_mz_tolerance > 0.0 && highest_intensity_peak_mz_tolerance <= 0.0)
        {
          LOG_INFO << "Corrected " << corrected_to_nearest_peak.size() << " precursor to a MS1 peak." << endl;
        }
        else if (highest_intensity_peak_mz_tolerance > 0.0)
        {
          LOG_INFO << "Corrected " << corrected_to_highest_intensity_peak.size() << " precursor to a MS1 peak." << endl;
        }
        else
        {
          LOG_WARN << "Output file 'out_csv': No data collected since 'nearest_peak:mz_tolerance' was not enabled. CSV will be empty." << endl;
        }
        writeHist(out_csv, deltaMZs, mzs, rts);
      }

      if (!in_feature.empty())
      {
        LOG_INFO << "Corrected " << corrected_to_nearest_feature.size() << " precursors to a feature." << endl;
      }

      return EXECUTION_OK;
    }

};

int main(int argc, const char ** argv)
{
  TOPPHiResPrecursorMassCorrector tool;
  return tool.main(argc, argv);
}

/// @endcond

