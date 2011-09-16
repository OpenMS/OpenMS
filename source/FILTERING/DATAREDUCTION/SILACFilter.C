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
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedIsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_randist.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <cmath>
#include <iostream>

using namespace std;

namespace OpenMS
{
  SILACFilter::SILACFilter(std::vector<DoubleReal> mass_separations, Int charge, DoubleReal model_deviation, Int isotopes_per_peptide)
  {
    mass_separations_ = mass_separations;     // mass shift(s) between peptides
    charge_ = charge;     // peptide charge
    model_deviation_ = model_deviation;   // allowed deviation from averegine model
    isotopes_per_peptide_ = isotopes_per_peptide;   // isotopic peaks per peptide

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

  SILACFilter::SILACFilter()
  {
  }

  SILACFilter::~SILACFilter()
  {

  }


  bool SILACFilter::extractMzShiftsAndIntensities(const MSSpectrum<Peak1D> &s, DoubleReal mz, DoubleReal picked_mz, const SILACFiltering &f)
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
        if (nearestPeakError < f.getPeakWidth(mz)/3)
        {
          deltaMZ = nearest_peak_mz - picked_mz;
        } else
        {
          deltaMZ = -1; // peak not found
        }


        // old
        /*
        DoubleReal expensive_MzShift = computeActualMzShift_(mz, expected_shift, getPeakWidth(mz));

        if (expensive_MzShift - deltaMZ > 0.01)
        {

         cout << "rt " << rt << " mz " << mz << " expected mz " << expected_shift << " deltaMZ " << deltaMZ << " expMZ " << expensive_MzShift << " diff " << expensive_MzShift - deltaMZ << endl;
        }
*/
        if ( deltaMZ < 0)
        {
          if (SILACFiltering::allow_missing_peaks_ == false)
          {
            return false;
          }
          else
          {
            // MISSING PEAK EXCEPTION
            // A missing intensity is allowed if (1) the user allowed it, (2) it's the last isotopic peak of a SILAC peptide and (3) it hasn't occured before.
            if (SILACFiltering::allow_missing_peaks_ == true && isotope == isotopes_per_peptide_ - 1 && missing_peak_seen_yet == false)
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
          exact_intensities_singlePeptide.push_back(gsl_spline_eval (SILACFiltering::spline_spl_, mz + deltaMZ, SILACFiltering::current_spl_));
        }
      }

      exact_shifts_.push_back(exact_shifts_singlePeptide);
      exact_intensities_.push_back(exact_intensities_singlePeptide);
      expected_shifts_.push_back(expected_shifts_singlePeptide);      // store expected_shifts for blacklisting
    }
    return true;
  }

  bool SILACFilter::extractMzShiftsAndIntensitiesPicked(const MSSpectrum<Peak1D> &s, DoubleReal mz, const SILACFiltering &f)
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
          cout << "exp_mz, mz, err, pw/2 " << next_peak_expected << " " << nearest_peak_mz << " " << nearestPeakError << " " << f.getPeakWidth(mz)/2.0 << endl;
        }

        // check if error is small enough
        if (nearestPeakError < f.getPeakWidth(mz)/2.0)
        {
          deltaMZ = nearest_peak_mz - mz;
        } else
        {
          deltaMZ = -1; // peak not found
        }

        if ( deltaMZ < 0)
        {
          if (SILACFiltering::allow_missing_peaks_ == false)
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
            if (SILACFiltering::allow_missing_peaks_ == true && isotope == isotopes_per_peptide_ - 1 && missing_peak_seen_yet == false)
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

  bool SILACFilter::extractMzShiftsAndIntensitiesPickedToPattern(const MSSpectrum<Peak1D> &s, DoubleReal mz, const SILACFiltering &f, SILACPattern &pattern)
  {
    if (!extractMzShiftsAndIntensitiesPicked(s, mz, f))
      return false;

    pattern.rt = s.getRT();
    pattern.mz = mz;
    pattern.charge = charge_;
    pattern.isotopes_per_peptide = (Int) isotopes_per_peptide_;
    pattern.intensities.insert(pattern.intensities.begin(), exact_intensities_.begin(), exact_intensities_.end());
    pattern.mass_shifts.insert(pattern.mass_shifts.begin(), mz_peptide_separations_.begin(), mz_peptide_separations_.end());
    return true;
  }

  bool SILACFilter::intensityFilter()
  {
    bool missing_peak_seen_yet = false;  // Did we encounter a missing peak in this SILAC pattern yet?
    for (Size peptide = 0; peptide != number_of_peptides_; ++peptide)
    {
      for (Size isotope = 0; isotope < isotopes_per_peptide_; ++isotope)
      {
        if (exact_intensities_[peptide][isotope] < SILACFiltering::intensity_cutoff_)
        {
          if (SILACFiltering::allow_missing_peaks_ == false)
          {
            return false;
          }
          else
          {
            // MISSING PEAK EXCEPTION
            // A missing intensity is allowed if (1) the user allowed it, (2) it's the last isotopic peak of a SILAC peptide and (3) it hasn't occured before.
            if (SILACFiltering::allow_missing_peaks_ == true && isotope == isotopes_per_peptide_ - 1 && missing_peak_seen_yet == false)
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

  bool SILACFilter::correlationFilter1(DoubleReal mz, const SILACFiltering &f)
  {
    bool missing_peak_seen_yet = false;

    for (Size peptide = 0; peptide != number_of_peptides_; ++peptide)
    {
      for (Size isotope2 = 1; isotope2 < isotopes_per_peptide_; ++isotope2)
      {
        std::vector<DoubleReal> intensities1;    // intensities in region around first peak of peptide
        std::vector<DoubleReal> intensities2;    // intensities in region around following peak
        DoubleReal mzWindow = 0.7 * f.getPeakWidth(mz);    // width of the window around m/z in which the correlation is calculated

        for (DoubleReal dmz = - mzWindow; dmz <= mzWindow; dmz += 0.2 * mzWindow)     // fill intensity vectors
        {
          DoubleReal intens1 = gsl_spline_eval(SILACFiltering::spline_spl_, mz + exact_shifts_[peptide][0] + dmz, SILACFiltering::current_spl_);
          DoubleReal intens2 = gsl_spline_eval(SILACFiltering::spline_spl_, mz + exact_shifts_[peptide][isotope2] + dmz, SILACFiltering::current_spl_);
          intensities1.push_back( intens1 );
          intensities2.push_back( intens2 );
        }

        DoubleReal intensityCorrelation = Math::pearsonCorrelationCoefficient( intensities1.begin(), intensities1.end(), intensities2.begin(), intensities2.end());    // calculate Pearson correlation coefficient

        if (intensityCorrelation < SILACFiltering::intensity_correlation_)
        {
          // MISSING PEAK EXCEPTION
          // A missing intensity is allowed if (1) the user allowed it, (2) one of the two peaks is the last isotopic peak of a SILAC peptide and (3) it hasn't occured before.
          if (SILACFiltering::allow_missing_peaks_ && (isotope2 == isotopes_per_peptide_ - 1) && (!missing_peak_seen_yet))
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

  bool SILACFilter::correlationFilter2(DoubleReal mz, const SILACFiltering &f)
  {
    if (!(number_of_peptides_ == 1 && exact_shifts_[0][0] == 0))     // If we are looking for single peptides, this filter is not needed.
    {
      for (Size peptide = 0; peptide < number_of_peptides_ - 1; ++peptide)
      {
        std::vector<DoubleReal> intensities3;    // intensities in region around monoisotopic peak
        std::vector<DoubleReal> intensities4;    // intensities in region around first peak of following peptide
        DoubleReal mzWindow = 0.7 * f. getPeakWidth(mz);    // width of the window around m/z in which the correlation is calculated

        for (DoubleReal dmz = - mzWindow; dmz <= mzWindow; dmz += 0.2 * mzWindow)     // fill intensity vectors
        {
          DoubleReal intens3 = gsl_spline_eval(SILACFiltering::spline_spl_, mz + exact_shifts_[0][0] + dmz, SILACFiltering::current_spl_);
          DoubleReal intens4 = gsl_spline_eval(SILACFiltering::spline_spl_, mz + exact_shifts_[peptide + 1][0] + dmz, SILACFiltering::current_spl_);
          intensities3.push_back( intens3 );
          intensities4.push_back( intens4 );
        }

        DoubleReal intensityCorrelation = Math::pearsonCorrelationCoefficient( intensities3.begin(), intensities3.end(), intensities4.begin(), intensities4.end());    // calculate Pearson correlation coefficient

        if (intensityCorrelation < SILACFiltering::intensity_correlation_)
        {
          return false;
        }
      }
    }
    return true;
  }

  bool SILACFilter::averageneFilter(DoubleReal mz)
  {
    bool missing_peak_seen_yet = false;

    // bool debug = abs(rt - 6653.3) < 0.1 && abs(mz - 668.83) < 0.01;
    bool debug = false;

    if (!IsotopeDistributionCache::getInstance()->isPrecalculated())
    {
      IsotopeDistributionCache::getInstance()->precalculate(20000.0, 1.0, 0.0, 0.0);
    }

    if (isotopes_per_peptide_ > 1)
    {
      for (Size peptide = 0; peptide != number_of_peptides_; ++peptide)
      {
        //IsotopeDistribution isoDistribution;    // isotope distribution of an averagene peptide
        //isoDistribution.estimateFromPeptideWeight((mz + exact_shifts_[peptide][0]) * charge_);    // mass of averagene peptide

        const TheoreticalIsotopePattern& pattern = IsotopeDistributionCache::getInstance()->getIsotopeDistribution((mz + exact_shifts_[peptide][0]) * charge_);

        DoubleReal averagineIntensity_mono = pattern.intensity[0];    // intensity of monoisotopic peak of the averagine model
        DoubleReal intensity_mono = exact_intensities_[peptide][0];    // intensity around the (potential) monoisotopic peak in the real data

        for (Size isotope = 1; isotope < isotopes_per_peptide_; ++isotope)
        {
          DoubleReal averagineIntensity = pattern.intensity[isotope];
          DoubleReal intensity = exact_intensities_[peptide][isotope];

          DoubleReal averageneRatio = averagineIntensity / averagineIntensity_mono;

          DoubleReal ratio = (intensity / intensity_mono) / averageneRatio;

          DoubleReal max_deviation = model_deviation_ / averageneRatio;

          if (debug)
          {
            cout << "Ratio: " << ratio << endl;
          }
          if (ratio > max_deviation || ratio < 1 / max_deviation) // Test for missing peak?
          {
            // MISSING PEAK EXCEPTION
            // A missing intensity is allowed if (1) the user allowed it, (2) one of the two peaks is the last isotopic peak of a SILAC peptide and (3) it hasn't occured before.
            if (SILACFiltering::allow_missing_peaks_ && (isotope == isotopes_per_peptide_ - 1) && (!missing_peak_seen_yet))
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

  bool SILACFilter::isSILACPattern_(const MSSpectrum<Peak1D> &s, DoubleReal mz, DoubleReal picked_mz, const SILACFiltering &f, MSSpectrum<Peak1D> &debug, SILACPattern &pattern)
  {
    current_mz_ = mz;

    Peak1D debug_peak;
    debug_peak.setMZ(mz);

    // EXACT m/z SHIFTS (Determine the actual shifts between peaks. Say 4 Th is the theoretic shift. In the experimental data it will be 4.0029 Th.)
    if (!extractMzShiftsAndIntensities(s, mz, picked_mz, f))
    {
      debug_peak.setIntensity(1);
      debug.push_back(debug_peak);
      return false;
    }

    // COMPLETE INTENSITY FILTER (Check that all of the intensities are above the cutoff.)
    if (!intensityFilter())
    {
      debug_peak.setIntensity(2);
      debug.push_back(debug_peak);
      return false;
    }

    // CORRELATION FILTER 1 (Check for every peptide that its mono-isotopic peak correlates with the following peaks)
    if (!correlationFilter1(mz, f))
    {
      debug_peak.setIntensity(3);
      debug.push_back(debug_peak);
      return false;
    }

    // CORRELATION FILTER 2 (Check that the monoisotopic peak of the light (unlabeled) peptide correlates with the mono-isotopic peak of the labeled peptides)
    if (!correlationFilter2(mz, f))
    {
      debug_peak.setIntensity(4);
      debug.push_back(debug_peak);
      return false;
    }

    // AVERAGINE FILTER (Check if realtive ratios confirm with an averagine model of all peptides.)
    if (!averageneFilter(mz))
    {
      debug_peak.setIntensity(5);
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
    if (!extractMzShiftsAndIntensitiesPicked(s, mz, f))
    {
      debug_peak.setIntensity(1);
      debug.push_back(debug_peak);
      return false;
    }

    // COMPLETE INTENSITY FILTER (Check that all of the intensities are above the cutoff.)
    if (!intensityFilter())
    {
      debug_peak.setIntensity(2);
      debug.push_back(debug_peak);
      return false;
    }

    // AVERAGINE FILTER (Check if realtive ratios confirm with an averagine model of all peptides.)
    if (!averageneFilter(mz))
    {
      debug_peak.setIntensity(3);
      debug.push_back(debug_peak);
      return false;
    }

    debug_peak.setIntensity(10);
    debug.push_back(debug_peak);

    return true;
  }



  DoubleReal SILACFilter::computeActualMzShift_(DoubleReal mz, DoubleReal expectedMzShift, DoubleReal maxMzDeviation)
  {
    if (expectedMzShift <= 0.0)
    {
      return 0;			// return 0 for potential monoisotopic peak
    }
    else
    {
      //--------------------------------------------------
      // compute autocorrelation
      //--------------------------------------------------

      std::vector<DoubleReal> akimaMz;
      DoubleReal stepwidth = maxMzDeviation / 30;

      // n must be a power of two for gsl fast fourier transformation; take next higher size for n, which is a power of two and fill the rest with zeros
      Size akimaMz_size = pow(2,(ceil(log((3*maxMzDeviation)/stepwidth)/log(2.0))));

      akimaMz.clear();
      akimaMz.resize(akimaMz_size, 0.0);

      // check to not leave experiment
      DoubleReal starting_offset = std::min(mz - SILACFiltering::mz_min_, maxMzDeviation);

      // calculate akima interpolation for region around mz (+- maxMzDeviation) and store in vector akimaMz
      // starting position: mz - maxMzDeviation
      // ending position: mz + maxMzDeviation
      // 19 steps, stepwidth = maxMzDeviation / 30 each
      // akima interpolation for position x
      Size i = 0;
      for (DoubleReal x = mz - starting_offset; x <= mz + maxMzDeviation; x += stepwidth)
      {
        akimaMz[i] = gsl_spline_eval_deriv2 (SILACFiltering::spline_aki_, x, SILACFiltering::current_aki_);
        ++i;
      }

      // calculate akima interpolation for region around mz + expectedMzShift (- 2 * maxMzDeviation, + maxMzDeviation) and store in vector akimaMzShift
      // starting position: mz - 2 * maxMzDeviation
      // ending position: mz + maxMzDeviation
      // 19 steps, maxMzDeviation / 30 each
      // akima interpolation for position (x + expectedMzShift)
      std::vector<DoubleReal> akimaMzShift(akimaMz_size, 0.0);
      i = 0;
      for (DoubleReal x = mz - starting_offset - maxMzDeviation; x <= mz + maxMzDeviation; x += stepwidth)
      {
        akimaMzShift[i] = gsl_spline_eval_deriv2 (SILACFiltering::spline_aki_, x + expectedMzShift, SILACFiltering::current_aki_);
        ++i;
      }

      // fourier transformation of the values
      gsl_fft_real_radix2_transform (&(*akimaMz.begin()), 1, akimaMz_size);
      gsl_fft_real_radix2_transform (&(*akimaMzShift.begin()), 1, akimaMz_size);

      // create vector "autoCorrelations" to store autocorrelations
      std::vector<DoubleReal> autoCorrelations;
      autoCorrelations.resize(akimaMz_size);

      // multiply the fourier transformed complex values with the complex conjugate
      // have a look at the GNU Scientific Library reference manual for a description of the data structure.
      autoCorrelations[0] = akimaMz[0] * akimaMzShift[0];     // special case for i = 0

      autoCorrelations[akimaMz_size / 2] = akimaMz[akimaMz_size / 2] * akimaMzShift[akimaMz_size / 2];      // special case for i = akimaMz_size / 2

      for (i = 1; i <= akimaMz_size / 2; ++i)
      {
        autoCorrelations[i] = akimaMz[i] * akimaMzShift[i] + akimaMz[akimaMz_size - i] * akimaMzShift[akimaMz_size - i];      // cases for (0 < i < akimaMz_size/ 2)
        autoCorrelations[akimaMz_size - i] = 0.0;     // for (akimaMz_size / 2 < i akimaMz_size) fill the vector "autoCorrelations" with zeros
      }

      // compute inverse fourier transformation
      gsl_fft_halfcomplex_radix2_inverse (&(*autoCorrelations.begin()), 1, akimaMz_size);

      autoCorrelations.resize(akimaMz_size / 2);     // cut the vector in half to erase zero entries


      //--------------------------------------------------
      // compute exact position
      //--------------------------------------------------

      // create two vectors for interpolation
      std::vector<DoubleReal> mzShifts(autoCorrelations.size(), 0.0);      // vector "mzShifts" contains potetial mz shifts (starting by expectedMzShift - maxMzDeviation, moving by maxMzDeviation / 18, ending by the end of vector "autocorrelations")
      std::vector<DoubleReal> correspondingAutoCorrelations(autoCorrelations.size(), 0.0);      // vector "correspondingAutoCorrelations" contains autocorrelations that corresponds to potetial mz shifht in vector "mzShifts"

      DoubleReal shift = expectedMzShift - maxMzDeviation;      // caculate first potential mz shift

      // fill vectors "mzShifts" and "correspondingAutoCorrelations"
      for (Size i = 0; i < autoCorrelations.size(); ++i)
      {
        mzShifts[i] = shift;
        correspondingAutoCorrelations[i] = autoCorrelations[i];
        shift += stepwidth;
      }

      // interpolate the current autocorrelation peak
      gsl_interp_accel* acc_correlation = gsl_interp_accel_alloc();
      gsl_spline* spline_correlation = gsl_spline_alloc(gsl_interp_cspline, mzShifts.size());
      gsl_spline_init(spline_correlation, &(*mzShifts.begin()), &(*correspondingAutoCorrelations.begin()), mzShifts.size());

      DoubleReal last_intensity_fwd = 0.0;
      DoubleReal current_intensity_fwd = 0.0;
      DoubleReal next_intensity_fwd = 0.0;
      DoubleReal last_intensity_bwd = 0.0;
      DoubleReal current_intensity_bwd = 0.0;
      DoubleReal next_intensity_bwd = 0.0;

      for (DoubleReal current_position = 0.0; current_position <= maxMzDeviation; current_position += stepwidth)
      {
        // search for first maximum in + direction and check preconditions
        if (current_position == 0.0)
        {
          last_intensity_fwd = gsl_spline_eval (spline_correlation, expectedMzShift + current_position - stepwidth, acc_correlation);
          current_intensity_fwd = gsl_spline_eval (spline_correlation, expectedMzShift + current_position, acc_correlation);
          next_intensity_fwd = gsl_spline_eval (spline_correlation, expectedMzShift + current_position + stepwidth, acc_correlation);
        }

        else
        {
          last_intensity_fwd = current_intensity_fwd;
          current_intensity_fwd = next_intensity_fwd;
          next_intensity_fwd = gsl_spline_eval (spline_correlation, expectedMzShift + current_position + stepwidth, acc_correlation);
        }

        // search for a current m/z shift larger than the expected one
        // coditions: intensity at position (mz + expectedMzShift + current_position) > intensity_cutoff (intensity calculated with akima interpolation based on "intensities_vec" from SILACFiltering)
        if (current_intensity_fwd > last_intensity_fwd && current_intensity_fwd > next_intensity_fwd && gsl_spline_eval (SILACFiltering::spline_aki_, mz + expectedMzShift + current_position, SILACFiltering::current_aki_) > SILACFiltering::intensity_cutoff_)
        {
          gsl_spline_free(spline_correlation);      // free interpolation object
          gsl_interp_accel_free(acc_correlation);     // free accelerator object
          return expectedMzShift + current_position;      // return exact position
        }

        // search for first maximum in - direction and check preconditions
        if (current_position == 0.0)
        {
          last_intensity_bwd = gsl_spline_eval(spline_correlation, expectedMzShift - current_position - stepwidth, acc_correlation);
          current_intensity_bwd = gsl_spline_eval(spline_correlation, expectedMzShift - current_position, acc_correlation);
          next_intensity_bwd = gsl_spline_eval(spline_correlation, expectedMzShift - current_position + stepwidth, acc_correlation);
        }

        else
        {
          next_intensity_bwd = current_intensity_bwd;
          current_intensity_bwd = last_intensity_bwd;
          last_intensity_bwd = gsl_spline_eval(spline_correlation, expectedMzShift - current_position - stepwidth, acc_correlation);
        }

        // search for an current m/z shift smaller than the expected one
        // conditions: intensity at position (mz + expectedMzShift - current_position) > intensity_cutoff (intensity calculated with akima interpolation based on "intensities_vec" from SILACFiltering)
        if (current_intensity_bwd > last_intensity_bwd && current_intensity_bwd > next_intensity_bwd && gsl_spline_eval (SILACFiltering::spline_aki_, mz + expectedMzShift - current_position, SILACFiltering::current_aki_) > SILACFiltering::intensity_cutoff_)
        {
          gsl_spline_free(spline_correlation);      // free interpolation object
          gsl_interp_accel_free(acc_correlation);     // free accelerator object
          return expectedMzShift - current_position;      // return exact position
        }
      }

      gsl_spline_free(spline_correlation);      // free interpolation object
      gsl_interp_accel_free(acc_correlation);     // free accelerator object
      return -1;      // return -1 if no autocorrelation exists for expectedMzShift
    }
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
