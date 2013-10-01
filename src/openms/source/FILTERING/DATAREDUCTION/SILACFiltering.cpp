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
// $Maintainer: Lars Nilse $
// $Authors: Steffen Sass, Holger Plattfaut $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/DATAREDUCTION/SILACFiltering.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFilter.h>
#include <OpenMS/MATH/MISC/LinearInterpolation.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>

#include <iostream>
#include <fstream>
#include <map>
#include <list>

using namespace std;

namespace OpenMS
{
  SILACFiltering::SpectrumInterpolation::SpectrumInterpolation(const MSSpectrum<>& s, const SILACFiltering& f)
  {
    vector<DoubleReal> mz, intensity;
    DoubleReal last_mz = s.begin()->getMZ();
    DoubleReal peak_width_cur = f.peak_width(last_mz);

    // Fill intensity and m/z vector for interpolation. Add zeros in the area with no data points to improve cubic spline fit
    for (MSSpectrum<>::ConstIterator mz_interpol_it = s.begin(); mz_interpol_it != s.end(); ++mz_interpol_it)
    {
      if (mz_interpol_it->getMZ() > last_mz + peak_width_cur) // If the mz gap is rather larger, fill in zeros. These addtional Stützstellen improve interpolation where no signal (i.e. data points) is.
      {
        for (DoubleReal current_mz = last_mz + peak_width_cur; current_mz < mz_interpol_it->getMZ() - peak_width_cur; current_mz += peak_width_cur)
        {
          mz.push_back(current_mz);
          intensity.push_back(0.0);
        }
      }
      if (mz_interpol_it->getMZ() > last_mz)
      {
        mz.push_back(mz_interpol_it->getMZ());
        intensity.push_back(mz_interpol_it->getIntensity());
      }
      last_mz = mz_interpol_it->getMZ();
    }

    // spline interpolation, used for exact ratio calculation (more accurate when real peak pairs are present)
    current_ = deprecated_gsl_interp_accel_alloc();
    spline_ = deprecated_gsl_spline_alloc(deprecated_wrapper_get_gsl_interp_cspline(), mz.size());
    deprecated_gsl_spline_init(spline_, &mz.front(), &intensity.front(), mz.size());
  }

  SILACFiltering::SpectrumInterpolation::~SpectrumInterpolation()
  {
    deprecated_gsl_interp_accel_free(current_);
    deprecated_gsl_spline_free(spline_);
  }

  SILACFiltering::SILACFiltering(MSExperiment<Peak1D>& exp, const PeakWidthEstimator::Result& peak_width, const DoubleReal intensity_cutoff, const String debug_filebase) :
    intensity_cutoff_(intensity_cutoff),
    exp_(exp),
    debug_filebase_(debug_filebase),
    peak_width(peak_width)
  {
  }

  void SILACFiltering::addFilter(SILACFilter& filter)
  {
    filters_.push_back(filter);
  }

  void SILACFiltering::pickSeeds_()
  {
    // perform peak picking
    PeakPickerHiRes picker;
    picker.setLogType(getLogType());
    Param param = picker.getParameters();
    param.setValue("ms1_only", DataValue("true"));
    param.setValue("signal_to_noise", 0.1);
    picker.setParameters(param);

    picker.pickExperiment(exp_, picked_exp_);

    if (debug_filebase_ != "")
    {
      MzMLFile mz_data_file;
      mz_data_file.store(debug_filebase_ + ".filtering.picked.mzML", picked_exp_);
    }

    // Initialize seeds map
    picked_exp_seeds_ = picked_exp_;
    for (Size i = 0; i != picked_exp_seeds_.size(); ++i)
    {
      picked_exp_seeds_[i].clear(false);
    }
  }

  void SILACFiltering::filterSeeds_()
  {
    startProgress(0, filters_.size(), "filtering seed data");

    UInt filter_id = 0;
    // Iterate over all filters
    for (vector<SILACFilter>::iterator filter_it = filters_.begin(); filter_it != filters_.end(); ++filter_it, ++filter_id)
    {
      setProgress(filter_it - filters_.begin());

      MSExperiment<Peak1D> exp_debug;

      UInt picked_rt_id = 0;
      // Iterate over all spectra of the experiment (iterate over rt)
      for (MSExperiment<Peak1D>::Iterator picked_rt_it = picked_exp_.begin(); picked_rt_it != picked_exp_.end(); ++picked_rt_it, ++picked_rt_id)
      {
        DoubleReal rt = picked_rt_it->getRT();

        MSSpectrum<Peak1D> debug;
        debug.setRT(rt);
        debug.setMSLevel(1);
        debug.setNativeID(String("debug-seed=") + picked_rt_id);

        // Iterate over the picked spectrum
        for (MSSpectrum<Peak1D>::Iterator picked_mz_it = picked_rt_it->begin(); picked_mz_it != picked_rt_it->end(); ++picked_mz_it) // iteration correct
        {
          DoubleReal picked_mz = picked_mz_it->getMZ();

          bool isSILAC = filter_it->isSILACPatternPicked_(*picked_rt_it, picked_mz, *this, debug);

          if (isSILAC)
          {
            Size spec_idx = picked_rt_it - picked_exp_.begin();
            picked_exp_seeds_[spec_idx].push_back(*picked_mz_it);
          }
        }

        exp_debug.addSpectrum(debug);
      }

      if (debug_filebase_ != "")
      {
        ChromatogramTools().convertSpectraToChromatograms(exp_debug, true);
        Int mass_separation = 0;
        if (filter_it->mass_separations_.size())
          mass_separation = filter_it->mass_separations_[0];
        MzMLFile().store(debug_filebase_ + ".filtering.seeds-filters:" +
                         filter_it->charge_ + ";" +
                         mass_separation + ";" +
                         filter_it->isotopes_per_peptide_ + ";" +
                         filter_it->model_deviation_ +
                         ".mzML", exp_debug);
      }
    }

    // sort complete spectrum so we can use range queries
    picked_exp_seeds_.sortSpectra(true);
    endProgress();

    if (debug_filebase_ != "")
    {
      MzMLFile mz_data_file;
      mz_data_file.store(debug_filebase_ + ".filtering.seeds.mzML", picked_exp_seeds_);
    }
  }

  void SILACFiltering::filterDataPoints()
  {
    pickSeeds_();
    filterSeeds_();

    startProgress(0, exp_.size(), "filtering raw data");

    UInt filter_id = 0;
    // Iterate over all filters
    for (vector<SILACFilter>::iterator filter_it = filters_.begin(); filter_it != filters_.end(); ++filter_it, ++filter_id)
    {
      MSExperiment<Peak1D> exp_debug;

      UInt rt_id = 0;
      // Iterate over all spectra of the experiment (iterate over rt)
      for (MSExperiment<Peak1D>::Iterator picked_seed_rt_it = picked_exp_seeds_.begin();
           picked_seed_rt_it != picked_exp_seeds_.end();
           ++picked_seed_rt_it, ++rt_id)
      {
        DoubleReal rt = picked_seed_rt_it->getRT(); // retention time of this spectrum

        MSExperiment<Peak1D>::Iterator picked_rt_it = picked_exp_.RTBegin(rt);
        MSExperiment<Peak1D>::Iterator rt_it = exp_.RTBegin(rt);

        // set progress
        // calculate with progress for the current rt run and progress for the filter run, each scaled by total numbers of filters
        setProgress((picked_seed_rt_it - picked_exp_seeds_.begin()) / filters_.size() + distance(filters_.begin(), filter_it) * picked_exp_seeds_.size() / filters_.size());

        MSSpectrum<Peak1D> debug;
        debug.setRT(rt);
        debug.setMSLevel(1);
        debug.setNativeID(String("debug-spline=") + rt_id);

        // spectra with less than 10 data points and less then two picked peaks are being ignored
        if (rt_it->size() >= 10 && picked_rt_it->size() > 1)
        {
          SpectrumInterpolation spec_inter(*rt_it, *this);

          // XXX: Workaround to catch duplicated peaks
          std::set<DoubleReal> seen_mz;

          // Iterate over the picked spectrum
          for (MSSpectrum<Peak1D>::Iterator picked_mz_it = picked_seed_rt_it->begin(); picked_mz_it != picked_seed_rt_it->end(); ++picked_mz_it) // iteration correct
          {
            DoubleReal picked_mz = picked_mz_it->getMZ();
            DoubleReal intensity = picked_mz_it->getIntensity();

            // XXX: Ignore duplicated peaks
            if (!seen_mz.insert(picked_mz).second)
              continue;

            //---------------------------------------------------------------
            // BLUNT INTENSITY FILTER (Just check that intensity at current m/z position is above the intensity cutoff)
            //---------------------------------------------------------------
            if (intensity < intensity_cutoff_)
            {
              continue;
            }

            // XXX: Extract peaks again
            SILACPattern pattern;
            if (!filter_it->extractMzShiftsAndIntensitiesPickedToPattern_(*picked_rt_it, picked_mz, *this, pattern))
              continue;

            DoubleReal peak_width_cur = peak_width(picked_mz);

            for (DoubleReal mz = picked_mz - peak_width_cur; mz < picked_mz + peak_width_cur; mz += 0.1 * peak_width_cur) // iteration correct
            {
              //--------------------------------------------------
              // BLACKLIST FILTER
              //--------------------------------------------------

              // iterate over the blacklist (Relevant blacklist entries are most likely among the last ones added.)
              multimap<DoubleReal, BlacklistEntry>::iterator blacklistStartCheck;
              multimap<DoubleReal, BlacklistEntry>::iterator blacklistEndCheck;

              if (blacklist.size() > 40) // Blacklist should be of certain size before we ckeck only parts of it.
              {
                blacklistStartCheck = blacklist.lower_bound(rt - 100);
                blacklistEndCheck = blacklist.lower_bound(rt);
              }
              else
              {
                blacklistStartCheck = blacklist.begin();
                blacklistEndCheck = blacklist.end();
              }

              bool isBlacklisted = false;

              for (multimap<DoubleReal, BlacklistEntry>::iterator blacklist_check_it = blacklistStartCheck; blacklist_check_it != blacklistEndCheck; ++blacklist_check_it)
              {
                Int charge = filter_it->getCharge();
                const vector<DoubleReal>& mass_separations = filter_it->getMassSeparations();

                // loop over the individual isotopic peaks of the SILAC pattern (and check if they are blacklisted)
                const vector<DoubleReal>& expectedMZshifts = filter_it->getExpectedMzShifts();

                for (vector<DoubleReal>::const_iterator expectedMZshifts_it = expectedMZshifts.begin(); expectedMZshifts_it != expectedMZshifts.end(); ++expectedMZshifts_it)
                {
                  bool inBlacklistEntry = blacklist_check_it->second.range.encloses(*expectedMZshifts_it + mz, rt);
                  bool exception = (charge == blacklist_check_it->second.charge)
                                   && (mass_separations == blacklist_check_it->second.mass_separations)
                                   && (fabs(*expectedMZshifts_it - blacklist_check_it->second.relative_peak_position) < 0.1);

                  if (inBlacklistEntry && !exception)
                  {
                    isBlacklisted = true;
                    break;
                  }
                }

                if (isBlacklisted)
                {
                  break;
                }
              }

              // Check the other filters only if current m/z and rt position is not blacklisted
              if (isBlacklisted == false)
              {
                if (filter_it->isSILACPattern_(*picked_rt_it, spec_inter, mz, picked_mz, *this, debug, pattern)) // Check if the mz at the given position is a SILAC pair
                {
                  //--------------------------------------------------
                  // FILLING THE BLACKLIST
                  //--------------------------------------------------

                  DoubleReal peak_width_cur = peak_width(mz);

                  // loop over the individual isotopic peaks of the SILAC pattern (and blacklist the area around them)
                  const vector<DoubleReal>& peak_positions = filter_it->getPeakPositions();

                  // Remember the charge and mass separations (since the blacklisting should not apply to filters of the same charge and mass separations).
                  Int charge = filter_it->getCharge();
                  const std::vector<DoubleReal>& mass_separations = filter_it->getMassSeparations();

                  for (vector<DoubleReal>::const_iterator peak_positions_it = peak_positions.begin(); peak_positions_it != peak_positions.end(); ++peak_positions_it)
                  {
                    DRange<2> blackArea; // area in the m/z-RT plane to be blacklisted
                    blackArea.setMinX(*peak_positions_it - 0.8 * peak_width_cur); // set min m/z position of area to be blacklisted
                    blackArea.setMaxX(*peak_positions_it + 0.8 * peak_width_cur); // set max m/z position of area to be blacklisted
                    blackArea.setMinY(rt - 10); // set min rt position of area to be blacklisted
                    blackArea.setMaxY(rt + 10); // set max rt position of area to be blacklisted

                    // Remember relative m/z shift (since the blacklisting should not apply to filters of the same points of the same relative m/z shift).
                    DoubleReal relative_peak_position = *peak_positions_it - mz;

                    // Does the new black area overlap with existing areas in the blacklist?
                    bool overlap = false;
                    // Does the current filter and relative peak position agree with the ones of the blacklist entry?
                    bool sameFilterAndPeakPosition = false;

                    multimap<DoubleReal, BlacklistEntry>::iterator blacklistStartFill;
                    multimap<DoubleReal, BlacklistEntry>::iterator blacklistEndFill;
                    if (blacklist.size() > 40) // Blacklist should be of certain size before we ckeck only parts of it.
                    {
                      blacklistStartFill = blacklist.lower_bound(rt - 100);
                      blacklistEndFill = blacklist.lower_bound(rt);
                    }
                    else
                    {
                      blacklistStartFill = blacklist.begin();
                      blacklistEndFill = blacklist.end();
                    }
                    for (multimap<DoubleReal, BlacklistEntry>::iterator blacklist_fill_it = blacklistStartFill; blacklist_fill_it != blacklistEndFill; ++blacklist_fill_it)
                    {
                      overlap = blackArea.isIntersected(blacklist_fill_it->second.range);
                      sameFilterAndPeakPosition = (charge == blacklist_fill_it->second.charge) && (mass_separations == blacklist_fill_it->second.mass_separations) && (abs(relative_peak_position - blacklist_fill_it->second.relative_peak_position) < 0.01);

                      if (overlap && sameFilterAndPeakPosition)
                      {
                        // If new and old entry intersect, simply update (or replace) the old one.
                        if (blackArea.minY() > (blacklist_fill_it->second.range).minY())
                        {
                          // no new min RT => no change of key necessary
                          (blacklist_fill_it->second.range).setMinX(min(blackArea.minX(), (blacklist_fill_it->second.range).minX()));
                          (blacklist_fill_it->second.range).setMaxX(max(blackArea.maxX(), (blacklist_fill_it->second.range).maxX()));
                          (blacklist_fill_it->second.range).setMaxY(max(blackArea.maxY(), (blacklist_fill_it->second.range).maxY()));
                        }
                        else
                        {
                          // new min RT => insert new BlacklistEntry and delete old one
                          DRange<2> mergedArea;
                          BlacklistEntry mergedEntry;
                          mergedArea.setMinX(min(blackArea.minX(), (blacklist_fill_it->second.range).minX()));
                          mergedArea.setMaxX(max(blackArea.maxX(), (blacklist_fill_it->second.range).maxX()));
                          mergedArea.setMinY(blackArea.minY());
                          mergedArea.setMaxY(max(blackArea.maxY(), (blacklist_fill_it->second.range).maxY()));
                          mergedEntry.range = mergedArea;
                          mergedEntry.charge = blacklist_fill_it->second.charge;
                          mergedEntry.mass_separations = blacklist_fill_it->second.mass_separations;
                          mergedEntry.relative_peak_position = blacklist_fill_it->second.relative_peak_position;

                          // Simply insert the new and erase the old map BlacklistEntry. We break out of the loop anyhow.
                          blacklist.insert(pair<DoubleReal, BlacklistEntry>(mergedEntry.range.minY(), mergedEntry));
                          blacklist.erase(blacklist_fill_it);
                        }

                        break;
                      }
                    }

                    if (!overlap)
                    {
                      // If new and none of the old entries intersect, add a new entry.
                      BlacklistEntry newEntry;
                      newEntry.range = blackArea;
                      newEntry.charge = charge;
                      newEntry.mass_separations = mass_separations;
                      newEntry.relative_peak_position = relative_peak_position;
                      blacklist.insert(pair<DoubleReal, BlacklistEntry>(newEntry.range.minY(), newEntry));
                    }
                  }

                  /*                  // DEBUG: save global blacklist as .csv
                  ofstream blacklistFile;
                  blacklistFile.open ("blacklist.csv");

                  for (map<DoubleReal,BlacklistEntry>::iterator blacklist_it = blacklist.begin(); blacklist_it != blacklist.end(); ++blacklist_it)
                  {
                    blacklistFile << rt << ", " << (blacklist_it->second.range).minX() << ", " << (blacklist_it->second.range).maxX() << ", " << (blacklist_it->second.range).minY() << ", " << (blacklist_it->second.range).maxY() << ", " << (blacklist_it->second.charge) << ", " << (blacklist_it->second.mass_separations[0]) << ", " << (blacklist_it->second.relative_peak_position) << endl;
                  }
                  blacklistFile.close();
*/
                }
              }
            }

            // XXX
            const UInt threshold_points = 4;
            if (pattern.points.size() > threshold_points)
              filter_it->elements_.push_back(pattern);
          }
        }

        exp_debug.addSpectrum(debug);
      }

      if (debug_filebase_ != "")
      {
        ChromatogramTools().convertSpectraToChromatograms(exp_debug, true);
        Int mass_separation = 0;
        if (filter_it->mass_separations_.size())
          mass_separation = filter_it->mass_separations_[0];
        MzMLFile().store(debug_filebase_ + ".filtering.spline-filters:" +
                         filter_it->charge_ + ";" +
                         mass_separation + ";" +
                         filter_it->isotopes_per_peptide_ + ";" +
                         filter_it->model_deviation_ +
                         ".mzML", exp_debug);
      }
    }

    endProgress();
  }

}
