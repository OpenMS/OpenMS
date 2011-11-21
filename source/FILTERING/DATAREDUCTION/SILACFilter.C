// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Lars Nilse, Steffen Sass, Holger Plattfaut $
// --------------------------------------------------------------------------


#include <OpenMS/FILTERING/DATAREDUCTION/SILACFilter.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFiltering.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedIsotopeFitter1D.h>
#include <OpenMS/FILTERING/DATAREDUCTION/IsotopeDistributionCache.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <cmath>
#include <iostream>

using namespace std;

namespace OpenMS
{
  SILACFilter::SILACFilter(std::vector<DoubleReal> mass_separations, Int charge, DoubleReal model_deviation, Int isotopes_per_peptide,
      DoubleReal intensity_cutoff, DoubleReal intensity_correlation, bool allow_missing_peaks)
    : mass_separations_(mass_separations),
      charge_(charge),
      model_deviation_(model_deviation),
      isotopes_per_peptide_(isotopes_per_peptide),
      intensity_cutoff_(intensity_cutoff),
      intensity_correlation_(intensity_correlation),
      allow_missing_peaks_(allow_missing_peaks),
      isotope_distribution_(20000.0, 1.0, 0.0, 0.0)
  {
    isotope_distance_ = 1.000495 / (DoubleReal)charge_;    // distance between isotopic peaks of a peptide [Th]
    number_of_peptides_ = (Int) mass_separations_.size() + 1;    // number of labelled peptides +1 [e.g. for SILAC triplet =3]
    
    // m/z shifts from mass shifts
    mz_peptide_separations_.push_back(0.0);
    for (std::vector<DoubleReal>::iterator it = mass_separations_.begin(); it != mass_separations_.end(); ++it)
    {
      mz_peptide_separations_.push_back(*it / (DoubleReal)charge_);
    }
    
    expected_mz_shifts_.clear();
    for (std::vector<DoubleReal>::iterator it = mz_peptide_separations_.begin(); it != mz_peptide_separations_.end(); ++it)
    {
      for (Size i = 0; i < isotopes_per_peptide_; i++)
      {
        expected_mz_shifts_.push_back(*it + i * isotope_distance_);
      }
    }
  }

  bool SILACFilter::extractMzShiftsAndIntensities_(const MSSpectrum<Peak1D> &s, const SILACFiltering::SpectrumInterpolation &si, DoubleReal mz, DoubleReal picked_mz, const SILACFiltering &f)
  {
    bool missing_peak_seen_yet = false;  // Did we encounter a missing peak in this SILAC pattern yet?

    exact_shifts_.clear();
    exact_intensities_.clear();
    expected_shifts_.clear();

    for (Size peptide = 0; peptide != number_of_peptides_; ++peptide) // loop over labelled peptides [e.g. for SILAC triplets: 0=light 1=medium 2=heavy]
    {
      std::vector<DoubleReal> exact_shifts_singlePeptide;
      std::vector<DoubleReal> exact_intensities_singlePeptide;
      std::vector<DoubleReal> expected_shifts_singlePeptide;

      for (Size isotope = 0; isotope < isotopes_per_peptide_; isotope++) // loop over isotopic peaks within a peptide [0=mono-isotopic peak etc.]
      {

        DoubleReal expected_shift = mz_peptide_separations_[peptide] + isotope * isotope_distance_;

        // calculate expected position of next peak
        DoubleReal next_peak_expected = picked_mz + expected_shift;

        // find peak (index) closest to expected position
        Size nearest_peak_idx = s.findNearest(next_peak_expected);

        // get actual position of closest peak
        DoubleReal nearest_peak_mz = s[nearest_peak_idx].getMZ();

        // calculate error between expected and actual position
        DoubleReal nearestPeakError = abs( nearest_peak_mz - next_peak_expected);

        DoubleReal deltaMZ;

        // check if error is small enough
        if (nearestPeakError < f.peak_width(mz) / 3)
        {
          deltaMZ = nearest_peak_mz - picked_mz;
        } else
        {
          deltaMZ = -1; // peak not found
        }


        if ( deltaMZ < 0)
        {
          if (allow_missing_peaks_ == false)
          {
            return false;
          }
          else
          {
            // MISSING PEAK EXCEPTION
            // A missing intensity is allowed if (1) the user allowed it, (2) it's the last isotopic peak of a SILAC peptide and (3) it hasn't occured before.
            if (allow_missing_peaks_ == true && isotope == isotopes_per_peptide_ - 1 && missing_peak_seen_yet == false)
            {
              missing_peak_seen_yet = true;
            }
            else
            {
              return false;
            }
          }
        }

        exact_shifts_singlePeptide.push_back( deltaMZ );

        expected_shifts_singlePeptide.push_back(mz_peptide_separations_[peptide] + isotope * isotope_distance_);      // store expected_shift for blacklisting

        if ( deltaMZ < 0 )
        {
          exact_intensities_singlePeptide.push_back( -1 );
        }
        else
        {
          exact_intensities_singlePeptide.push_back(si(mz + deltaMZ));
        }
      }

      exact_shifts_.push_back(exact_shifts_singlePeptide);
      exact_intensities_.push_back(exact_intensities_singlePeptide);
      expected_shifts_.push_back(expected_shifts_singlePeptide);      // store expected_shifts for blacklisting
    }
    return true;
  }

  bool SILACFilter::extractMzShiftsAndIntensitiesPicked_(const MSSpectrum<Peak1D> &s, DoubleReal mz, const SILACFiltering &f)
  {
    //bool debug = abs(rt - 6653.3) < 0.1 && abs(mz - 668.83) < 0.01;
    bool debug = false;

    bool missing_peak_seen_yet = false;  // Did we encounter a missing peak in this SILAC pattern yet?

    exact_shifts_.clear();
    exact_intensities_.clear();
    expected_shifts_.clear();

    if (debug)
    {
      cout << "n-peptides: " << number_of_peptides_ << endl;
    }

    for (Size peptide = 0; peptide != number_of_peptides_; ++peptide) // loop over labelled peptides [e.g. for SILAC triplets: 0=light 1=medium 2=heavy]
    {
      std::vector<DoubleReal> exact_shifts_singlePeptide;
      std::vector<DoubleReal> exact_intensities_singlePeptide;
      std::vector<DoubleReal> expected_shifts_singlePeptide;

      for (Size isotope = 0; isotope < isotopes_per_peptide_; isotope++) // loop over isotopic peaks within a peptide [0=mono-isotopic peak etc.]
      {

        DoubleReal expected_shift = mz_peptide_separations_[peptide] + isotope * isotope_distance_;

        // calculate expected position of next peak
        DoubleReal next_peak_expected = mz + expected_shift;

        // find peak (index) closest to expected position
        Size nearest_peak_idx = s.findNearest(next_peak_expected);

        // get actual position and intensity of closest peak
        DoubleReal nearest_peak_mz = s[nearest_peak_idx].getMZ();
        DoubleReal nearest_peak_intensity = s[nearest_peak_idx].getIntensity();


        // calculate error between expected and actual position
        DoubleReal nearestPeakError = abs( nearest_peak_mz - next_peak_expected);

        DoubleReal deltaMZ;

        if(debug)
        {
          cout << "exp_mz, mz, err, pw/2 " << next_peak_expected << " " << nearest_peak_mz << " " << nearestPeakError << " " << f.peak_width(mz) / 2.0 << endl;
        }

        // check if error is small enough
        if (nearestPeakError < f.peak_width(mz) / 2.0)
        {
          deltaMZ = nearest_peak_mz - mz;
        } else
        {
          deltaMZ = -1; // peak not found
        }

        if ( deltaMZ < 0)
        {
          if (allow_missing_peaks_ == false)
          {
            if(debug)
            {
              cout << "Missing Peak!" << endl;
            }
            return false;
          }
          else
          {
            // MISSING PEAK EXCEPTION
            // A missing intensity is allowed if (1) the user allowed it, (2) it's the last isotopic peak of a SILAC peptide and (3) it hasn't occured before.
            if (allow_missing_peaks_ == true && isotope == isotopes_per_peptide_ - 1 && missing_peak_seen_yet == false)
            {
              missing_peak_seen_yet = true;
            }
            else
            {
              return false;
            }
          }
        }

        exact_shifts_singlePeptide.push_back( deltaMZ );

        expected_shifts_singlePeptide.push_back(mz_peptide_separations_[peptide] + isotope * isotope_distance_);      // store expected_shift for blacklisting

        if ( deltaMZ < 0 )
        {
          exact_intensities_singlePeptide.push_back( -1 );
        }
        else
        {
          exact_intensities_singlePeptide.push_back(nearest_peak_intensity);
        }
      }

      exact_shifts_.push_back(exact_shifts_singlePeptide);
      exact_intensities_.push_back(exact_intensities_singlePeptide);
      expected_shifts_.push_back(expected_shifts_singlePeptide);      // store expected_shifts for blacklisting
    }

    if(debug)
    {
      cout << "Success!" << endl;
    }
    return true;
  }

  bool SILACFilter::extractMzShiftsAndIntensitiesPickedToPattern_(const MSSpectrum<Peak1D> &s, DoubleReal mz, const SILACFiltering &f, SILACPattern &pattern)
  {
    if (!extractMzShiftsAndIntensitiesPicked_(s, mz, f))
      return false;

    pattern.rt = s.getRT();
    pattern.mz = mz;
    pattern.charge = charge_;
    pattern.isotopes_per_peptide = (Int) isotopes_per_peptide_;
    pattern.intensities.insert(pattern.intensities.begin(), exact_intensities_.begin(), exact_intensities_.end());
    pattern.mass_shifts.insert(pattern.mass_shifts.begin(), mz_peptide_separations_.begin(), mz_peptide_separations_.end());
    return true;
  }

  bool SILACFilter::intensityFilter_()
  {
    bool missing_peak_seen_yet = false;  // Did we encounter a missing peak in this SILAC pattern yet?
    for (Size peptide = 0; peptide != number_of_peptides_; ++peptide)
    {
      for (Size isotope = 0; isotope < isotopes_per_peptide_; ++isotope)
      {
        if (exact_intensities_[peptide][isotope] < intensity_cutoff_)
        {
          if (allow_missing_peaks_ == false)
          {
            return false;
          }
          else
          {
            // MISSING PEAK EXCEPTION
            // A missing intensity is allowed if (1) the user allowed it, (2) it's the last isotopic peak of a SILAC peptide and (3) it hasn't occured before.
            if (allow_missing_peaks_ == true && isotope == isotopes_per_peptide_ - 1 && missing_peak_seen_yet == false)
            {
              missing_peak_seen_yet = true;
            }
            else
            {
              return false;
            }
          }
        }
      }
    }
    return true;
  }

  bool SILACFilter::correlationFilter1_(const SILACFiltering::SpectrumInterpolation &si, DoubleReal mz, const SILACFiltering &f)
  {
    bool missing_peak_seen_yet = false;

    for (Size peptide = 0; peptide != number_of_peptides_; ++peptide)
    {
      for (Size isotope2 = 1; isotope2 < isotopes_per_peptide_; ++isotope2)
      {
        std::vector<DoubleReal> intensities1;    // intensities in region around first peak of peptide
        std::vector<DoubleReal> intensities2;    // intensities in region around following peak
        DoubleReal mzWindow = 0.7 * f.peak_width(mz);    // width of the window around m/z in which the correlation is calculated

        for (DoubleReal dmz = - mzWindow; dmz <= mzWindow; dmz += 0.2 * mzWindow)     // fill intensity vectors
        {
          DoubleReal intens1 = si(mz + exact_shifts_[peptide][0] + dmz);
          DoubleReal intens2 = si(mz + exact_shifts_[peptide][isotope2] + dmz);
          intensities1.push_back( intens1 );
          intensities2.push_back( intens2 );
        }

        DoubleReal intensityCorrelation = Math::pearsonCorrelationCoefficient( intensities1.begin(), intensities1.end(), intensities2.begin(), intensities2.end());    // calculate Pearson correlation coefficient

        if (intensityCorrelation < intensity_correlation_)
        {
          // MISSING PEAK EXCEPTION
          // A missing intensity is allowed if (1) the user allowed it, (2) one of the two peaks is the last isotopic peak of a SILAC peptide and (3) it hasn't occured before.
          if (allow_missing_peaks_ && (isotope2 == isotopes_per_peptide_ - 1) && (!missing_peak_seen_yet))
          {
            missing_peak_seen_yet = true;
          }
          else
          {
            return false;

          }
        }
      }
    }
    return true;
  }

  bool SILACFilter::correlationFilter2_(const SILACFiltering::SpectrumInterpolation &si, DoubleReal mz, const SILACFiltering &f)
  {
    if (!(number_of_peptides_ == 1 && exact_shifts_[0][0] == 0))     // If we are looking for single peptides, this filter is not needed.
    {
      for (Size peptide = 0; peptide < number_of_peptides_ - 1; ++peptide)
      {
        std::vector<DoubleReal> intensities3;    // intensities in region around monoisotopic peak
        std::vector<DoubleReal> intensities4;    // intensities in region around first peak of following peptide
        DoubleReal mzWindow = 0.7 * f.peak_width(mz);    // width of the window around m/z in which the correlation is calculated

        for (DoubleReal dmz = - mzWindow; dmz <= mzWindow; dmz += 0.2 * mzWindow)     // fill intensity vectors
        {
          DoubleReal intens3 = si(mz + exact_shifts_[0][0] + dmz);
          DoubleReal intens4 = si(mz + exact_shifts_[peptide + 1][0] + dmz);
          intensities3.push_back( intens3 );
          intensities4.push_back( intens4 );
        }

        DoubleReal intensityCorrelation = Math::pearsonCorrelationCoefficient( intensities3.begin(), intensities3.end(), intensities4.begin(), intensities4.end());    // calculate Pearson correlation coefficient

        if (intensityCorrelation < intensity_correlation_)
        {
          return false;
        }
      }
    }
    return true;
  }

  bool SILACFilter::averageneFilter_(DoubleReal mz)
  {
    bool missing_peak_seen_yet = false;

    bool debug = false;
    //debug = abs(mz - 1788) < 1;

    if (isotopes_per_peptide_ > 1)
    {
      for (Size peptide = 0; peptide != number_of_peptides_; ++peptide)
      {
        //IsotopeDistribution isoDistribution;    // isotope distribution of an averagene peptide
        //isoDistribution.estimateFromPeptideWeight((mz + exact_shifts_[peptide][0]) * charge_);    // mass of averagene peptide

        const TheoreticalIsotopePattern& pattern = isotope_distribution_.getIsotopeDistribution((mz + exact_shifts_[peptide][0]) * charge_);

        DoubleReal averagine_mono = pattern.intensity[0];    // intensity of monoisotopic peak of the averagine model
        DoubleReal intensity_mono = exact_intensities_[peptide][0];    // intensity around the (potential) monoisotopic peak in the real data
        DoubleReal ratio_mono = intensity_mono / averagine_mono;

        for (Size isotope = 1; isotope < isotopes_per_peptide_; ++isotope)
        {
          DoubleReal averagine_current = pattern.intensity[isotope];
          DoubleReal intensity_current = exact_intensities_[peptide][isotope];
          DoubleReal ratio_current = intensity_current / averagine_current;

          DoubleReal ratio = ratio_current / ratio_mono;

          if (debug)
          {
            cout << "Ratio: " << ratio << endl;
          }
          if (ratio > model_deviation_ || ratio < 1 / model_deviation_) // Test for missing peak?
          {
            // MISSING PEAK EXCEPTION
            // A missing intensity is allowed if (1) the user allowed it, (2) one of the two peaks is the last isotopic peak of a SILAC peptide and (3) it hasn't occured before.
            if (allow_missing_peaks_ && (isotope == isotopes_per_peptide_ - 1) && (!missing_peak_seen_yet))
            {
              missing_peak_seen_yet = true;
            }
            else
            {
              //cout << "Missing Peak in averagine filter!" << endl;
              return false;
            }
          }
        }
      }
    }
    return true;
  }

  bool SILACFilter::isSILACPattern_(const MSSpectrum<Peak1D> &s, const SILACFiltering::SpectrumInterpolation &si, DoubleReal mz, DoubleReal picked_mz, const SILACFiltering &f, MSSpectrum<Peak1D> &debug, SILACPattern &pattern)
  {
    current_mz_ = mz;

    Peak1D debug_peak;
    debug_peak.setMZ(mz);

    // EXACT m/z SHIFTS (Determine the actual shifts between peaks. Say 4 Th is the theoretic shift. In the experimental data it will be 4.0029 Th.)
    if (!extractMzShiftsAndIntensities_(s, si, mz, picked_mz, f))
    {
      debug_peak.setIntensity(1);
      debug.push_back(debug_peak);
      return false;
    }

    // COMPLETE INTENSITY FILTER (Check that all of the intensities are above the cutoff.)
    if (!intensityFilter_())
    {
      debug_peak.setIntensity(2);
      debug.push_back(debug_peak);
      return false;
    }

    // CORRELATION FILTER 1 (Check for every peptide that its mono-isotopic peak correlates with the following peaks)
    if (!correlationFilter1_(si, mz, f))
    {
      debug_peak.setIntensity(3);
      debug.push_back(debug_peak);
      return false;
    }

    // CORRELATION FILTER 2 (Check that the monoisotopic peak of the light (unlabeled) peptide correlates with the mono-isotopic peak of the labeled peptides)
    if (!correlationFilter2_(si, mz, f))
    {
      debug_peak.setIntensity(4);
      debug.push_back(debug_peak);
      return false;
    }

    // ALL FILTERS PASSED => CREATE DATAPOINT
    SILACPoint newElement;    // Raw data point at this particular RT and m/z passed all filters. Store it for further clustering.
    newElement.rt = s.getRT();
    newElement.mz = mz;
    pattern.points.push_back(newElement);

    debug_peak.setIntensity(10);
    debug.push_back(debug_peak);

    return true;
  }

  bool SILACFilter::isSILACPatternPicked_(const MSSpectrum<Peak1D> &s, DoubleReal mz, const SILACFiltering &f, MSSpectrum<Peak1D> &debug)
  {
    current_mz_ = mz;

    Peak1D debug_peak;
    debug_peak.setMZ(mz);

    // EXACT m/z SHIFTS (Determine the actual shifts between peaks. Say 4 Th is the theoretic shift. In the experimental data it will be 4.0029 Th.)
    if (!extractMzShiftsAndIntensitiesPicked_(s, mz, f))
    {
      debug_peak.setIntensity(1);
      debug.push_back(debug_peak);
      return false;
    }

    // COMPLETE INTENSITY FILTER (Check that all of the intensities are above the cutoff.)
    if (!intensityFilter_())
    {
      debug_peak.setIntensity(2);
      debug.push_back(debug_peak);
      return false;
    }

    // AVERAGINE FILTER (Check if realtive ratios confirm with an averagine model of all peptides.)
    if (!averageneFilter_(mz))
    {
      debug_peak.setIntensity(3);
      debug.push_back(debug_peak);
      return false;
    }

    debug_peak.setIntensity(10);
    debug.push_back(debug_peak);

    return true;
  }

  std::vector<DoubleReal> SILACFilter::getPeakPositions()
	{
    peak_positions_.clear();
    for (Size peptide = 0; peptide != number_of_peptides_; ++peptide)
		{
      for (Size isotope = 0; isotope < isotopes_per_peptide_; ++isotope)
			{
        peak_positions_.push_back(current_mz_ + expected_shifts_[peptide][isotope]);
			}
		}
    return peak_positions_;
  }
	
  const std::vector<DoubleReal>& SILACFilter::getExpectedMzShifts()
	{
    return expected_mz_shifts_;
  }

  std::vector<SILACPattern>& SILACFilter::getElements()
  {
    return elements_;
  }

  Int SILACFilter::getCharge()
  {
    return charge_;
  }

  std::vector<DoubleReal>& SILACFilter::getMassSeparations()
  {
    return mass_separations_;
  }
}
