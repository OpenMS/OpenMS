// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
//

#include <OpenMS/PROCESSING/CALIBRATION/PrecursorCorrection.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <iomanip>
#include <fstream>

using namespace std;
using namespace OpenMS;

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

    void PrecursorCorrection::writeHist(const String& out_csv,
                                        const vector<double> & delta_mzs,
                                        const vector<double> & mzs,
                                        const vector<double> & rts)
    {
      //cout << "writing data" << endl;
      ofstream csv_file(out_csv.c_str());
      csv_file << setprecision(9);

      // header
      csv_file << ListUtils::concatenate(ListUtils::create<String>(PrecursorCorrection::csv_header), "\t") << "\n";

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

      for (Size scan = 0; scan != exp.size(); ++scan)
      {
        // skip non-tandem mass spectra
        if (exp[scan].getMSLevel() != 2 || exp[scan].getPrecursors().empty()) continue;

        // extract precursor / MS2 information
        const double pc_mz = exp[scan].getPrecursors()[0].getMZ();
        const double rt = exp[scan].getRT();
        const int pc_charge = exp[scan].getPrecursors()[0].getCharge();

        for (Size f = 0; f != features.size(); ++f)
        {
          // feature  is incompatible if believe_charge is set and charges don't match
          if (believe_charge && features[f].getCharge() != pc_charge)
          {
            continue;
          }
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

    bool PrecursorCorrection::overlaps_(const Feature& feature,
                                        const double rt,
                                        const double pc_mz,
                                        const double rt_tolerance)
    {
      if (feature.getConvexHulls().empty())
      {
        OPENMS_LOG_WARN << "HighResPrecursorMassCorrector warning: at least one feature has no convex hull - omitting feature for matching" << std::endl;
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

