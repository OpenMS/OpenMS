// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Steffen Sass, Holger Plattfaut $
// --------------------------------------------------------------------------
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFiltering.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFilter.h>
#include <OpenMS/MATH/MISC/LinearInterpolation.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <iostream>
#include <fstream>

using namespace std;

namespace OpenMS
{
  DoubleReal SILACFiltering::mz_stepwidth = 0;
  DoubleReal SILACFiltering::intensity_cutoff = 0;
  DoubleReal SILACFiltering::intensity_correlation = 0;
  bool SILACFiltering::allow_missing_peaks = true;
  gsl_interp_accel* SILACFiltering::current_aki = 0;
  gsl_interp_accel* SILACFiltering::current_spl = 0;
  gsl_spline* SILACFiltering::spline_aki = 0;
  gsl_spline* SILACFiltering::spline_spl = 0;
  Int SILACFiltering::feature_id = 0;
  DoubleReal SILACFiltering::mz_min = 0;

  SILACFiltering::SILACFiltering(MSExperiment<Peak1D>& exp_, DoubleReal mz_stepwidth_, DoubleReal intensity_cutoff_, DoubleReal intensity_correlation_, bool allow_missing_peaks_) : exp(exp_)
  {
    mz_stepwidth = mz_stepwidth_;
    intensity_cutoff = intensity_cutoff_;
    intensity_correlation = intensity_correlation_;
    allow_missing_peaks = allow_missing_peaks_;
  }

  void SILACFiltering::addFilter(SILACFilter& filter)
  {
    filters.push_back(&filter);
  }

  SILACFiltering::~SILACFiltering()
  {
  }

  void SILACFiltering::filterDataPoints()
  {
    startProgress(0, exp.size(), "filtering raw data");

    vector<DataPoint> data;

    mz_min = exp.getMinMZ();      // get lowest m/z value

    // Iterate over all filters
    for (list<SILACFilter*>::iterator filter_it = filters.begin(); filter_it != filters.end(); ++filter_it)
    {
      // Iterate over all spectra of the experiment (iterate over rt)
      for (MSExperiment<Peak1D>::Iterator rt_it = exp.begin(); rt_it != exp.end(); ++rt_it)
      {
        // set progress
        // calculate with progress for the current rt run and progress for the filter run, each scaled by total numbers of filters
        setProgress((rt_it - exp.begin()) / filters.size() + distance(filters.begin(), filter_it) * exp.size() / filters.size());

        Size number_data_points = rt_it->size();    // number of (m/z, intensity) data points in this spectrum

        DoubleReal rt = rt_it->getRT();    // retention time of this spectrum

        // spectra with less than 10 data points are being ignored
        if (number_data_points >= 10)
        {
          // filter MS1 spectra
          // read one spectrum into GSL structure
          vector<DoubleReal> mz_vec;
          vector<DoubleReal> intensity_vec;
          mz_min = rt_it->begin()->getMZ();
          DoubleReal last_mz = rt_it->begin()->getMZ();

          // INTERPOLATION (Akima and Spline interpolation in order to have intensities at any m/z.)
          // Fill intensity and m/z vector for interpolation. Add zeros in the area with no data points to improve cubic spline fit
          for (MSSpectrum<>::Iterator mz_it = rt_it->begin(); mz_it != rt_it->end(); ++mz_it)
          {
            if (mz_it->getMZ() > last_mz + 2 * mz_stepwidth) // If the mz gap is rather larger, fill in zeros. These addtional Stützstellen improve interpolation where no signal (i.e. data points) is.
            {
              for (DoubleReal current_mz = last_mz + 2 * mz_stepwidth; current_mz < mz_it->getMZ() - 2 * mz_stepwidth; current_mz += mz_stepwidth)
              {
                mz_vec.push_back(current_mz);
                intensity_vec.push_back(0.0);
              }
            }
            mz_vec.push_back(mz_it->getMZ());
            intensity_vec.push_back(mz_it->getIntensity());
            last_mz = mz_it->getMZ();
          }

          // akima interpolation, returns 0 in regions with no raw data points
          current_aki = gsl_interp_accel_alloc();
          spline_aki = gsl_spline_alloc(gsl_interp_akima, mz_vec.size());
          gsl_spline_init(spline_aki, &*mz_vec.begin(), &*intensity_vec.begin(), mz_vec.size());

          // spline interpolation, used for exact ratio calculation (more accurate when real peak pairs are present)
          current_spl = gsl_interp_accel_alloc();
          spline_spl = gsl_spline_alloc(gsl_interp_cspline, mz_vec.size());
          gsl_spline_init(spline_spl, &*mz_vec.begin(), &*intensity_vec.begin(), mz_vec.size());

          MSSpectrum<>::Iterator mz_it = rt_it->begin();

          last_mz = mz_it->getMZ();
          ++mz_it;

          // Iterate over the spectrum with a step width that is oriented on the raw data point positions (iterate over mz)
          for ( ; mz_it != rt_it->end(); ++mz_it) // iteration correct
          {
            // We do not move with mz_stepwidth over the spline fit, but with about a third of the local mz differences
            for (DoubleReal mz = last_mz; mz < mz_it->getMZ(); mz += (abs(mz_it->getMZ() - last_mz)) / 3)
            {
              //---------------------------------------------------------------
              // BLUNT INTENSITY FILTER (Just check that intensity at current m/z position is above the intensity cutoff)
              //---------------------------------------------------------------

              if (gsl_spline_eval (spline_aki, mz, current_aki) < intensity_cutoff)
              {
                continue;
              }


              //--------------------------------------------------
              // BLACKLIST FILTER
              //--------------------------------------------------

              bool isBlacklisted = false;

              // iterate over the blacklist (Relevant blacklist entries are most likely among the last ones added.)
              multimap<DoubleReal, BlacklistEntry>::iterator blacklistStart;
              multimap<DoubleReal, BlacklistEntry>::iterator blacklistEnd;
              if (blacklist.size() > 40)    // Blacklist should be of certain size before we run ckeck only parts of it.
              {
                blacklistStart = blacklist.lower_bound(rt - 100);
                blacklistEnd = blacklist.lower_bound(rt);
              }
              else
              {
                blacklistStart = blacklist.begin();
                blacklistEnd = blacklist.end();
              }
              for (multimap<DoubleReal, BlacklistEntry>::iterator blacklist_it = blacklistStart; blacklist_it != blacklistEnd; ++blacklist_it)
              {
                
                Int charge = (*filter_it)->getCharge();
                vector<DoubleReal> mass_separations = (*filter_it)->getMassSeparations();
                
                // loop over the individual isotopic peaks of the SILAC pattern (and check if they are blacklisted)
                const vector<DoubleReal>& expectedMZshifts = (*filter_it)->getExpectedMZshifts();
                for (vector<DoubleReal>::const_iterator expectedMZshifts_it = expectedMZshifts.begin(); expectedMZshifts_it != expectedMZshifts.end(); ++expectedMZshifts_it)
                {
                  bool inBlacklistEntry = blacklist_it->second.range.encloses(*expectedMZshifts_it + mz, rt);
                  bool exception = (charge == blacklist_it->second.charge) && (mass_separations == blacklist_it->second.mass_separations) && (abs(*expectedMZshifts_it - blacklist_it->second.relative_peak_position)<0.01);
                  
                  if (inBlacklistEntry && !exception )
                  {
                    isBlacklisted = true;
                    break;
                  }
                }
              }
              

              // Check the other filters only if current m/z and rt position is not blacklisted
              if (isBlacklisted == false)
              {
                if ((*filter_it)->isSILACPattern(rt, mz))      // Check if the mz at the given position is a SILAC pair
                {
                  //--------------------------------------------------
                  // FILLING THE BLACKLIST
                  //--------------------------------------------------
                  
                  DoubleReal peak_width = SILACFilter::getPeakWidth(mz);
                  
                  // loop over the individual isotopic peaks of the SILAC pattern (and blacklist the area around them)
                  const vector<DoubleReal>& peak_positions = (*filter_it)->getPeakPositions();
                  for (vector<DoubleReal>::const_iterator peak_positions_it = peak_positions.begin(); peak_positions_it != peak_positions.end(); ++peak_positions_it)
                  {
                    DRange<2> blackArea;    // area in the m/z-RT plane to be blacklisted
                    blackArea.setMinX(*peak_positions_it - 0.8 * peak_width);     // set min m/z position of area to be blacklisted
                    blackArea.setMaxX(*peak_positions_it + 0.8 * peak_width);     // set max m/z position of area to be blacklisted
                    blackArea.setMinY(rt - 10);     // set min rt position of area to be blacklisted
                    blackArea.setMaxY(rt + 10);     // set max rt position of area to be blacklisted
                    
                    // Remember the charge, mass separations and relative m/z shift (since the blacklisting should not apply to filters of the same charge, mass separations and points of the same relative m/z shift).
                    Int charge = (*filter_it)->charge;
                    std::vector<DoubleReal> mass_separations;
                    mass_separations.insert(mass_separations.begin(), (*filter_it)->mass_separations.begin(), (*filter_it)->mass_separations.end());
                    DoubleReal relative_peak_position = *peak_positions_it - mz; // or mz_it->getMZ() ??
                    
                    // Does the new black area overlap with existing areas in the blacklist?
                    bool overlap = false;
                    // Does the current filter and relative peak position agree with the ones of the blacklist entry?
                    bool sameFilterAndPeakPosition = false;
                    
                    multimap<DoubleReal, BlacklistEntry>::iterator blacklistStart;
                    multimap<DoubleReal, BlacklistEntry>::iterator blacklistEnd;
                    if (blacklist.size() > 40)    // Blacklist should be of certain size before we run ckeck only parts of it.
                        {
                          blacklistStart = blacklist.lower_bound(rt - 100);
                          blacklistEnd = blacklist.lower_bound(rt);
                        }
                        else
                        {
                          blacklistStart = blacklist.begin();
                          blacklistEnd = blacklist.end();
                        }
                    for (multimap<DoubleReal, BlacklistEntry>::iterator blacklist_it = blacklistStart; blacklist_it != blacklistEnd; ++blacklist_it)
                    {
                      overlap = blackArea.isIntersected(blacklist_it->second.range);
                      sameFilterAndPeakPosition = (charge == blacklist_it->second.charge) && (mass_separations == blacklist_it->second.mass_separations) && (abs(relative_peak_position - blacklist_it->second.relative_peak_position)<0.01);
                      
                      if (overlap && sameFilterAndPeakPosition)
                      {
                        // If new and old entry intersect, simply update (or replace) the old one.
                        if (blackArea.minY() > (blacklist_it->second.range).minY())
                         {
                           // no new min RT => no change of key necessary
                           (blacklist_it->second.range).setMinX(min(blackArea.minX(),(blacklist_it->second.range).minX()));
                           (blacklist_it->second.range).setMaxX(max(blackArea.maxX(),(blacklist_it->second.range).maxX()));
                           (blacklist_it->second.range).setMaxY(max(blackArea.maxY(),(blacklist_it->second.range).maxY()));
                         }
                         else
                         {
                           // new min RT => insert new BlacklistEntry and delete old one
                           DRange<2> mergedArea;
                           BlacklistEntry mergedEntry;
                           mergedArea.setMinX(min(blackArea.minX(), (blacklist_it->second.range).minX()));
                           mergedArea.setMaxX(max(blackArea.maxX(), (blacklist_it->second.range).maxX()));
                           mergedArea.setMinY(blackArea.minY());
                           mergedArea.setMaxY(max(blackArea.maxY(), (blacklist_it->second.range).maxY()));
                           mergedEntry.range = mergedArea;
                           mergedEntry.charge = blacklist_it->second.charge;
                           mergedEntry.mass_separations = blacklist_it->second.mass_separations;
                           mergedEntry.relative_peak_position = blacklist_it->second.relative_peak_position;
                           
                           // Simply insert the new and erase the old map BlacklistEntry. We break out of the loop anyhow.
                           blacklist.insert(pair<DoubleReal, BlacklistEntry>(mergedEntry.range.minY(), mergedEntry));
                           blacklist.erase(blacklist_it);
                         }
                        
                        break;
                      }
                    }
                    
                    if ( !overlap )
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
                  
/*                  // DEBUG: save global blacklist
                  ofstream blacklistFile;
                  blacklistFile.open ("blacklist.csv");
                  for (map<DoubleReal,BlacklistEntry>::iterator blacklist_it = blacklist.begin(); blacklist_it != blacklist.end(); ++blacklist_it)
                  {
                    blacklistFile << rt << ", " << (blacklist_it->second.range).minX() << ", " << (blacklist_it->second.range).maxX() << ", " << (blacklist_it->second.range).minY() << ", " << (blacklist_it->second.range).maxY() << ", " << (blacklist_it->second.charge) << ", " << (blacklist_it->second.mass_separations[0]) << ", " << (blacklist_it->second.relative_peak_position) << endl;
                  }
                  blacklistFile.close();
*/

                  ++feature_id;
                }
              }
            }
            
            last_mz = mz_it->getMZ();
          }
      }

        // Clear the interpolations
        gsl_spline_free(spline_aki);
        gsl_interp_accel_free(current_aki);
        gsl_spline_free(spline_spl);
        gsl_interp_accel_free(current_spl);
      }
    }

    endProgress();
  }

}
