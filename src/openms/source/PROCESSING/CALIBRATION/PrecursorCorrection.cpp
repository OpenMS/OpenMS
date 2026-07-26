// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
//

#include <OpenMS/PROCESSING/CALIBRATION/PrecursorCorrection.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <utility>
#include <vector>

using namespace std;
using namespace OpenMS;

namespace
{
  /**
    @brief The region in which a feature accepts a precursor.

    This is the bounding box of the feature's convex hull, extended by @p rt_tolerance in
    retention time and by 0.01 Th in m/z. All four edges are inclusive.
  */
  struct FeatureBox
  {
    double rt_lo, rt_hi, mz_lo, mz_hi;

    /// Spelled as the negation of "outside", so that a NaN coordinate is accepted exactly as
    /// DBoundingBox::encloses() accepts it (neither comparison fires) rather than rejected.
    bool contains(double rt, double mz) const
    {
      return !(rt < rt_lo || rt > rt_hi || mz < mz_lo || mz > mz_hi);
    }

    /// Whether this box can be placed on an ordered m/z sweep axis at all
    bool sortable() const
    {
      return std::isfinite(mz_lo) && std::isfinite(mz_hi);
    }
  };

  FeatureBox makeFeatureBox(const OpenMS::Feature& feature, double rt_tolerance)
  {
    // Deliberately built with the same DBoundingBox operations the per-pair check always used,
    // so the region is identical for every input - including the normalisation setMin()/setMax()
    // perform when a tolerance is negative, or when the feature has no convex hull at all.
    OpenMS::DBoundingBox<2> bb = feature.getConvexHull().getBoundingBox();
    const OpenMS::DPosition<2> extend_rt(rt_tolerance, 0.01);
    bb.setMin(bb.minPosition() - extend_rt);
    bb.setMax(bb.maxPosition() + extend_rt);
    return FeatureBox{bb.minPosition()[0], bb.maxPosition()[0],
                      bb.minPosition()[1], bb.maxPosition()[1]};
  }
} // namespace

namespace OpenMS
{

   const std::string PrecursorCorrection::csv_header = "RT,uncorrectedMZ,correctedMZ,deltaMZ";


   void PrecursorCorrection::getPrecursors(const MSExperiment & exp,
                                           vector<Precursor> & precursors,
                                           vector<double> & precursors_rt,
                                           vector<Size> & precursor_scan_index)
    {
      for (Size i = 0; i != exp.size(); ++i)
      {
        const vector<Precursor>& pcs = exp[i].getPrecursors();
        if (pcs.empty()) { continue; }
        vector<double> pcs_rt(pcs.size(), exp[i].getRT());
        copy(pcs.begin(), pcs.end(), back_inserter(precursors));
        copy(pcs_rt.begin(), pcs_rt.end(), back_inserter(precursors_rt));
        precursor_scan_index.push_back(i);
      }
    }

    void PrecursorCorrection::writeHist(const std::string& out_csv,
                                        const vector<double> & delta_mzs,
                                        const vector<double> & mzs,
                                        const vector<double> & rts)
    {
      //cout << "writing data" << endl;
      ofstream csv_file(out_csv.c_str());
      csv_file << setprecision(9);

      // header
      csv_file << ListUtils::concatenate(ListUtils::create<std::string>(PrecursorCorrection::csv_header), "\t") << "\n";

      // entries
      for (vector<double>::const_iterator it = delta_mzs.begin(); it != delta_mzs.end(); ++it)
      {
        UInt index = it - delta_mzs.begin();
        csv_file << rts[index] << "\t" << mzs[index] << "\t" << mzs[index] + *it  << "\t" << *it << "\n";
      }
      csv_file.close();
    }

     set<Size> PrecursorCorrection::correctToNearestMS1Peak(MSExperiment & exp,
                                                            double mz_tolerance,
                                                            bool ppm,
                                                            vector<double> & delta_mzs,
                                                            vector<double> & mzs,
                                                            vector<double> & rts)
    {
      set<Size> corrected_precursors;
      // load experiment and extract precursors
      vector<Precursor> precursors;  // precursor
      vector<double> precursors_rt;  // RT of precursor MS2 spectrum
      vector<Size> precursor_scan_index;
      getPrecursors(exp, precursors, precursors_rt, precursor_scan_index);

      for (Size i = 0; i != precursors_rt.size(); ++i)
      {
        // get precursor rt
        double rt = precursors_rt[i];

        // get precursor MZ
        double mz = precursors[i].getMZ();

        //cout << rt << " " << mz << endl;

        // get precursor spectrum
        MSExperiment::ConstIterator rt_it = exp.RTBegin(rt - 1e-8);

        // store index of MS2 spectrum
        UInt precursor_spectrum_idx = rt_it - exp.begin();

        // get parent (MS1) of precursor spectrum
        rt_it = exp.getPrecursorSpectrum(rt_it);

        if (rt_it == exp.end() 
        || rt_it->getMSLevel() != 1)
        {
          OPENMS_LOG_WARN << "Warning: no MS1 spectrum for this precursor" << endl;
          continue;          
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
            OPENMS_LOG_WARN << "Error: index is referencing different precursors in original and picked spectrum." << endl;
          }

          // cout << mz << " -> " << nearest_peak_mz << endl;
          double delta_mz = nearest_peak_mz - mz;
          delta_mzs.push_back(delta_mz);
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
    set<Size> PrecursorCorrection::correctToHighestIntensityMS1Peak(MSExperiment & exp,
                                                                    double mz_tolerance,
                                                                    bool ppm,
                                                                    vector<double> & delta_mzs,
                                                                    vector<double> & mzs,
                                                                    vector<double> & rts)
    {
      set<Size> corrected_precursors;
      // load experiment and extract precursors
      vector<Precursor> precursors;  // precursor
      vector<double> precursors_rt;  // RT of precursor MS2 spectrum
      vector<Size> precursor_scan_index;
      getPrecursors(exp, precursors, precursors_rt, precursor_scan_index);
      int count_error_highest_intenstiy = 0;

      for (Size i = 0; i != precursors_rt.size(); ++i)
      {
        double rt = precursors_rt[i]; // get precursor rt        
        double mz = precursors[i].getMZ(); // get precursor MZ

        // retrieves iterator of the MS2 fragment spectrum
        MSExperiment::ConstIterator rt_it = exp.RTBegin(rt - 1e-8);

        // store index of MS2 spectrum
        UInt precursor_spectrum_idx = rt_it - exp.begin();

        // get parent (MS1) of precursor spectrum
        rt_it = exp.getPrecursorSpectrum(rt_it);

        if (rt_it == exp.end() 
        || rt_it->getMSLevel() != 1)
        {
          OPENMS_LOG_WARN << "Warning: no MS1 spectrum for this precursor" << endl;
          continue;
        }

        // get tolerance window and index of highest peak
        std::pair<double,double> tolerance_window = Math::getTolWindow(mz, mz_tolerance, ppm);
        int highest_peak_idx = rt_it->findHighestInWindow(mz, mz-tolerance_window.first, tolerance_window.second-mz);

        // no MS1 precursor peak in +- tolerance window found
        if (highest_peak_idx == -1)
        {
          count_error_highest_intenstiy += 1;
          continue;
        }

        // get actual position and intensity of highest intensity peak
        double highest_peak_mz = (*rt_it)[highest_peak_idx].getMZ();
        double highest_peak_int = (*rt_it)[highest_peak_idx].getIntensity();

        // cout << mz << " -> " << nearest_peak_mz << endl;
        double delta_mz = highest_peak_mz - mz;
        delta_mzs.push_back(delta_mz);
        mzs.push_back(mz);
        rts.push_back(rt);
        // correct entries
        Precursor corrected_prec = precursors[i];
        corrected_prec.setMZ(highest_peak_mz);
        corrected_prec.setIntensity(highest_peak_int);
        exp[precursor_spectrum_idx].getPrecursors()[0] = corrected_prec;
        corrected_precursors.insert(precursor_spectrum_idx);
      }

      if (count_error_highest_intenstiy != 0)
      {
        OPENMS_LOG_INFO << "Correction to the highest intensity peak failed " 
           << count_error_highest_intenstiy 
           << " times because of missing peaks in the MS1. No changes were applied in these cases." 
           << std::endl;
      }

      return corrected_precursors;
    }

    set<Size> PrecursorCorrection::correctToNearestFeature(const FeatureMap& features,
                                                           MSExperiment & exp,
                                                           double rt_tolerance_s,
                                                           double mz_tolerance,
                                                           bool ppm,
                                                           bool believe_charge,
                                                           bool keep_original,
                                                           bool all_matching_features,
                                                           int max_trace,
                                                           int debug_level)
    {
      set<Size> corrected_precursors;

      // for each precursor/MS2 find all features that are in the given tolerance window (bounding box + rt tolerances)
      // if believe_charge is set, only add features that match the precursor charge
      map<Size, set<Size> > scan_idx_to_feature_idx;

      // MS2 scan indices ordered by precursor m/z. The experiment itself is never reordered.
      // The m/z is carried alongside the index so the sort compares plain doubles in contiguous
      // memory rather than dereferencing spectrum -> precursor list -> m/z on every comparison.
      // A non-finite precursor m/z has no place on the sweep axis - NaN would additionally
      // violate the strict weak ordering std::sort requires - so those spectra are kept aside
      // and compared exhaustively below.
      vector<pair<double, Size> > scan_by_mz;
      vector<Size> unsorted_scans;
      scan_by_mz.reserve(exp.size()); // upper bound; avoids reallocating while collecting
      for (Size scan = 0; scan != exp.size(); ++scan)
      {
        // skip non-tandem mass spectra
        if (exp[scan].getMSLevel() != 2 || exp[scan].getPrecursors().empty()) continue;
        const double pc_mz = exp[scan].getPrecursors()[0].getMZ();
        if (std::isfinite(pc_mz))
        {
          scan_by_mz.push_back(make_pair(pc_mz, scan));
        }
        else
        {
          unsorted_scans.push_back(scan);
        }
      }
      sort(scan_by_mz.begin(), scan_by_mz.end(),
           [](const pair<double, Size>& a, const pair<double, Size>& b) { return a.first < b.first; });

      // Pre-compute every feature's match region once. ConvexHull2D::getBoundingBox() is not
      // cached, so evaluating it per (spectrum, feature) pair used to dominate the runtime.
      // Skipped entirely when there is nothing to match against.
      vector<FeatureBox> boxes;
      vector<pair<double, Size> > feature_by_mz_lo; // ordered by the start of the m/z interval
      vector<Size> unsorted_features;               // non-finite m/z interval: compared directly
      if (!scan_by_mz.empty() || !unsorted_scans.empty())
      {
        boxes.reserve(features.size());
        feature_by_mz_lo.reserve(features.size());
        Size features_without_hull(0);
        for (Size f = 0; f != features.size(); ++f)
        {
          if (features[f].getConvexHulls().empty()) { ++features_without_hull; }
          boxes.push_back(makeFeatureBox(features[f], rt_tolerance_s));
          if (boxes[f].sortable())
          {
            feature_by_mz_lo.push_back(make_pair(boxes[f].mz_lo, f));
          }
          else
          {
            unsorted_features.push_back(f);
          }
        }

        // Reported once here instead of once per (spectrum, feature) pair, as the per-pair check did.
        // Such a feature still takes part in the sweep: its box normalises to a degenerate interval
        // at the numeric maximum, which no real precursor m/z can reach. That matches what the
        // exhaustive implementation did, warning wording included.
        if (features_without_hull != 0)
        {
          OPENMS_LOG_WARN << "HighResPrecursorMassCorrector warning: " << features_without_hull
                          << " feature(s) have no convex hull - omitting them for matching" << endl;
        }

        sort(feature_by_mz_lo.begin(), feature_by_mz_lo.end(),
             [](const pair<double, Size>& a, const pair<double, Size>& b) { return a.first < b.first; });
      }

      // Sweep a line over increasing precursor m/z: a feature becomes active when its m/z
      // interval starts and is retired once it ends, so every precursor is only compared
      // against the features whose m/z interval actually contains it.
      size_t overlap_checks(0);
      vector<Size> active;
      Size next_feature(0);
      for (Size i = 0; i != scan_by_mz.size(); ++i)
      {
        const double pc_mz = scan_by_mz[i].first;
        const Size scan = scan_by_mz[i].second;

        // extract precursor / MS2 information
        const double rt = exp[scan].getRT();
        const int pc_charge = exp[scan].getPrecursors()[0].getCharge();

        // activate intervals starting at or before pc_mz (must happen before retiring below)
        while (next_feature != feature_by_mz_lo.size() && feature_by_mz_lo[next_feature].first <= pc_mz)
        {
          active.push_back(feature_by_mz_lo[next_feature].second);
          ++next_feature;
        }

        Size kept(0);
        for (Size a = 0; a != active.size(); ++a)
        {
          const Size f = active[a];

          // interval ended; pc_mz only increases from here, so retire the feature for good
          if (boxes[f].mz_hi < pc_mz) { continue; }
          active[kept++] = f; // stays active regardless of the charge filter below

          // feature  is incompatible if believe_charge is set and charges don't match
          if (believe_charge && features[f].getCharge() != pc_charge)
          {
            continue;
          }
          // check if precursor/MS2 position overlap with feature
          if (boxes[f].contains(rt, pc_mz))
          {
            scan_idx_to_feature_idx[scan].insert(f);
          }
          ++overlap_checks;
        }
        active.resize(kept);

        // features that cannot be placed on the sweep axis (empty for well-formed input)
        for (Size u = 0; u != unsorted_features.size(); ++u)
        {
          const Size f = unsorted_features[u];
          if (believe_charge && features[f].getCharge() != pc_charge)
          {
            continue;
          }
          if (boxes[f].contains(rt, pc_mz))
          {
            scan_idx_to_feature_idx[scan].insert(f);
          }
          ++overlap_checks;
        }
      }

      // precursors that cannot be placed on the sweep axis are compared against every feature,
      // exactly as the exhaustive implementation did (empty for well-formed input)
      for (Size i = 0; i != unsorted_scans.size(); ++i)
      {
        const Size scan = unsorted_scans[i];
        const double pc_mz = exp[scan].getPrecursors()[0].getMZ();
        const double rt = exp[scan].getRT();
        const int pc_charge = exp[scan].getPrecursors()[0].getCharge();

        for (Size f = 0; f != features.size(); ++f)
        {
          if (believe_charge && features[f].getCharge() != pc_charge)
          {
            continue;
          }
          if (boxes[f].contains(rt, pc_mz))
          {
            scan_idx_to_feature_idx[scan].insert(f);
          }
          ++overlap_checks;
        }
      }

      if (debug_level > 0)
      {
        OPENMS_LOG_INFO << "Number of candidate overlap checks: " << overlap_checks
                        << " (comparing every precursor against every feature would be "
                        << (static_cast<UInt64>(scan_by_mz.size()) + static_cast<UInt64>(unsorted_scans.size()))
                           * static_cast<UInt64>(features.size())
                        << ")" << endl;
        OPENMS_LOG_INFO << "Number of precursors with overlapping features: " << scan_idx_to_feature_idx.size() << endl;
      }

      // filter sets to retain compatible features:
      // if precursor_mz = feature_mz + n * feature_charge (+/- mz_tolerance) a feature is compatible, others are removed from the set
      for (map<Size, set<Size> >::iterator it = scan_idx_to_feature_idx.begin(); it != scan_idx_to_feature_idx.end(); ++it)
      {
        const Size scan = it->first;
        const double pc_mz = exp[scan].getPrecursors()[0].getMZ();
        const double mz_tolerance_da = ppm ? pc_mz * mz_tolerance * 1e-6  : mz_tolerance;

        // Note: This is the "delete while iterating" pattern
        for (set<Size>::iterator sit = it->second.begin(); sit != it->second.end(); )
        {
          if (!compatible_(features[*sit], pc_mz, mz_tolerance_da, max_trace))
          {
            sit = it->second.erase(sit);
          }
          else
          {
            ++sit;
          }
        }
      }

      // remove entries with no compatible features (empty sets).
      // Note: This is the "delete while iterating" pattern
      for (map<Size, set<Size> >::iterator it = scan_idx_to_feature_idx.begin(); it != scan_idx_to_feature_idx.end(); )
      {
        if (it->second.empty())
        {
          it = scan_idx_to_feature_idx.erase(it);
        }
        else
        {
          ++it;
        }
      }

      if (debug_level > 0)
      {
        OPENMS_LOG_INFO << "Number of precursors with compatible features: " << scan_idx_to_feature_idx.size() << endl;
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
          // Note: This is the "delete while iterating" pattern
          for (set<Size>::iterator sit = it->second.begin(); sit != it->second.end(); )
          {
            if (sit != best_feature)
            {
              sit = it->second.erase(sit);
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

    bool PrecursorCorrection::overlaps_(const Feature& feature,
                                        const double rt,
                                        const double pc_mz,
                                        const double rt_tolerance)
    {
      if (feature.getConvexHulls().empty())
      {
        OPENMS_LOG_WARN << "HighResPrecursorMassCorrector warning: at least one feature has no convex hull - omitting feature for matching" << std::endl;
      }

      // bounding box extended by the retention time tolerance; correctToNearestFeature()
      // evaluates the very same region, but pre-computed once per feature
      return makeFeatureBox(feature, rt_tolerance).contains(rt, pc_mz);
    }

    bool PrecursorCorrection::compatible_(const Feature& feature,
                                          double pc_mz,
                                          double mz_tolerance,
                                          Size max_trace_number,
                                          int debug_level)
    {
      const int f_charge = feature.getCharge();
      const double f_mz = feature.getMZ();
      double trace = Math::round((pc_mz - f_mz) / (Constants::C13C12_MASSDIFF_U / f_charge)); // isotopic trace number at precursor mz
      double mass_error = fabs(pc_mz - (f_mz + trace * (Constants::C13C12_MASSDIFF_U / f_charge)));

      if (mass_error < mz_tolerance && (trace < max_trace_number + 0.01))
      {
        if (debug_level > 1)
        {
          OPENMS_LOG_INFO << "trace: " << (int)(trace + 0.5) << " feature_rt:" << feature.getRT() << " feature_mz:" << feature.getMZ() << " precursor_mz:" << pc_mz << endl;
        }
        return true;
      }
      else
      {
        return false;
      }
    }
}

