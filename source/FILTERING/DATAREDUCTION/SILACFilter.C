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
#include <OpenMS/DATASTRUCTURES/DataPoint.h>
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
#include <cmath>

namespace OpenMS
{
  SILACFilter::SILACFilter(std::vector<DoubleReal> mass_separations_, Int charge_, DoubleReal model_deviation_, Int isotopes_per_peptide_)
  {
    mass_separations = mass_separations_;     // mass shift(s) between peptides
    charge = charge_;     // peptide charge
    model_deviation = model_deviation_;   // allowed deviation from averegine model
    isotopes_per_peptide = isotopes_per_peptide_;   // isotopic peaks per peptide

    isotope_distance = 1.000495 / (DoubleReal)charge;    // distance between isotopic peaks of a peptide [Th]
    numberOfPeptides = mass_separations.size();    // number of labelled peptides +1 [e.g. for SILAC triplet =3]
    
    // m/z shifts from mass shifts
    mz_peptide_separations.push_back(0.0);
    for (std::vector<DoubleReal>::iterator it = mass_separations.begin(); it != mass_separations.end(); ++it)
    {
      mz_peptide_separations.push_back(*it / (DoubleReal)charge);
    }
    
    expectedMZshifts.clear();
    for (std::vector<DoubleReal>::iterator it = mz_peptide_separations.begin(); it != mz_peptide_separations.end(); ++it)
    {
      for (Int i=0; i<isotopes_per_peptide; i++)
      {
        expectedMZshifts.push_back(*it + i*isotope_distance);
      }
    }
    
    
    
  }

  SILACFilter::~SILACFilter()
  {

  }

  bool SILACFilter::isSILACPattern(DoubleReal rt, DoubleReal mz)
  {
    current_mz = mz;
    bool missing_peak_seen_yet = false;  // Did we encounter a missing peak in this SILAC pattern yet?


    //---------------------------------------------------------------
    // EXACT m/z SHIFTS (Determine the actual shifts between peaks. Say 4 Th is the theoretic shift. In the experimental data it will be 4.0029 Th.)
    //---------------------------------------------------------------
    exact_shifts.clear();
    exact_intensities.clear();
    expected_shifts.clear();

    for (Int peptide = 0; peptide <= numberOfPeptides; peptide++) // loop over labelled peptides [e.g. for SILAC triplets: 0=light 1=medium 2=heavy]
    {
      std::vector<DoubleReal> exact_shifts_singlePeptide;
      std::vector<DoubleReal> exact_intensities_singlePeptide;
      std::vector<DoubleReal> expected_shifts_singlePeptide;

      for (Int isotope = 0; isotope < isotopes_per_peptide; isotope++) // loop over isotopic peaks within a peptide [0=mono-isotopic peak etc.]
      {
        DoubleReal deltaMZ = computeActualMzShift(mz, mz_peptide_separations[peptide] + isotope * isotope_distance, getPeakWidth(mz));
        exact_shifts_singlePeptide.push_back( deltaMZ );

        expected_shifts_singlePeptide.push_back(mz_peptide_separations[peptide] + isotope * isotope_distance);      // store expected_shift for blacklisting

        if ( deltaMZ < 0 )
        {
          exact_intensities_singlePeptide.push_back( -1 );
        }
        else
        {
          exact_intensities_singlePeptide.push_back( gsl_spline_eval (SILACFiltering::spline_spl, mz + deltaMZ, SILACFiltering::current_spl) );
        }
      }

      exact_shifts.push_back(exact_shifts_singlePeptide);
      exact_intensities.push_back(exact_intensities_singlePeptide);
      expected_shifts.push_back(exact_shifts_singlePeptide);      // store expected_shifts for blacklisting
    }


    //---------------------------------------------------------------
    // COMPLETE INTENSITY FILTER (Check that all of the intensities are above the cutoff.)
    //---------------------------------------------------------------
    missing_peak_seen_yet = false;  // Did we encounter a missing peak in this SILAC pattern yet?
    for (Int peptide = 0; peptide <= numberOfPeptides; peptide++)
    {
      for (Int isotope = 0; isotope < isotopes_per_peptide; isotope++)
      {
        if (exact_intensities[peptide][isotope] < SILACFiltering::intensity_cutoff)
        {
          if (SILACFiltering::allow_missing_peaks == false)
          {
            return false;
          }
          else
          {
            // MISSING PEAK EXCEPTION
            // A missing intensity is allowed if (1) the user allowed it, (2) it's the last isotopic peak of a SILAC peptide and (3) it hasn't occured before.
            if (SILACFiltering::allow_missing_peaks == true && isotope == isotopes_per_peptide - 1 && missing_peak_seen_yet == false)
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


    //---------------------------------------------------------------
    // CORRELATION FILTER 1 (Check for every peptide that peak one correlates to following peaks of the same peptide)
    //---------------------------------------------------------------
    missing_peak_seen_yet = false;
    for (Int peptide = 0; peptide <= numberOfPeptides; peptide++)
    {
      for (Int isotope2 = 1; isotope2 < isotopes_per_peptide; isotope2++)
      {
        std::vector<DoubleReal> intensities1;    // intensities in region around first peak of peptide
        std::vector<DoubleReal> intensities2;    // intensities in region around following peak
        DoubleReal mzWindow = 0.7 * getPeakWidth(mz);    // width of the window around m/z in which the correlation is calculated

        for (DoubleReal dmz = - mzWindow; dmz <= mzWindow; dmz += 0.2 * mzWindow)     // fill intensity vectors
        {
          DoubleReal intens1 = gsl_spline_eval(SILACFiltering::spline_spl, mz + exact_shifts[peptide][0] + dmz, SILACFiltering::current_spl);
          DoubleReal intens2 = gsl_spline_eval(SILACFiltering::spline_spl, mz + exact_shifts[peptide][isotope2] + dmz, SILACFiltering::current_spl);
          intensities1.push_back( intens1 );
          intensities2.push_back( intens2 );
        }

        DoubleReal intensityCorrelation = Math::pearsonCorrelationCoefficient( intensities1.begin(), intensities1.end(), intensities2.begin(), intensities2.end());    // calculate Pearson correlation coefficient
        if ( intensityCorrelation < SILACFiltering::intensity_correlation )
        {
          // MISSING PEAK EXCEPTION
          // A missing intensity is allowed if (1) the user allowed it, (2) one of the two peaks is the last isotopic peak of a SILAC peptide and (3) it hasn't occured before.
          if (SILACFiltering::allow_missing_peaks && (isotope2 == isotopes_per_peptide - 1) && (!missing_peak_seen_yet))
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


    //---------------------------------------------------------------
    // CORRELATION FILTER 2 (Check that the monoisotopic peak correlates to every first peak of following peptides)
    //---------------------------------------------------------------
    for (Int peptide = 0; peptide < numberOfPeptides; peptide++)
    {
      std::vector<DoubleReal> intensities3;    // intensities in region around monoisotopic peak
      std::vector<DoubleReal> intensities4;    // intensities in region around first peak of following peptide
      DoubleReal mzWindow = 0.7 * getPeakWidth(mz);    // width of the window around m/z in which the correlation is calculated

      for (DoubleReal dmz = - mzWindow; dmz <= mzWindow; dmz += 0.2 * mzWindow)     // fill intensity vectors
      {
        DoubleReal intens3 = gsl_spline_eval(SILACFiltering::spline_spl, mz + exact_shifts[0][0] + dmz, SILACFiltering::current_spl);
        DoubleReal intens4 = gsl_spline_eval(SILACFiltering::spline_spl, mz + exact_shifts[peptide+1][0] + dmz, SILACFiltering::current_spl);
        intensities3.push_back( intens3 );
        intensities4.push_back( intens4 );
      }

      DoubleReal intensityCorrelation = Math::pearsonCorrelationCoefficient( intensities3.begin(), intensities3.end(), intensities4.begin(), intensities4.end());    // calculate Pearson correlation coefficient
      if ( intensityCorrelation < SILACFiltering::intensity_correlation )
      {
        return false;
      }
    }


    //---------------------------------------------------------------
    // AVERAGINE FILTER (Check if realtive ratios confirm with an averagine model of all peptides.)
    //---------------------------------------------------------------
    missing_peak_seen_yet = false;
    if (isotopes_per_peptide > 1)
    {
      for (Int peptide = 0; peptide <= numberOfPeptides; peptide++)
      {
        IsotopeDistribution isoDistribution;    // isotope distribution of an averagine peptide
        isoDistribution.estimateFromPeptideWeight((mz + exact_shifts[peptide][0]) * charge);    // mass of averagine peptide
        DoubleReal averagineIntensity_mono = isoDistribution.getContainer()[0].second;    // intensity of monoisotopic peak of the averagine model
        DoubleReal intensity_mono = exact_intensities[peptide][0];    // intensity around the (potential) monoisotopic peak in the real data

        for (Int isotope = 1; isotope < isotopes_per_peptide; isotope++)
        {
          DoubleReal averagineIntensity = isoDistribution.getContainer()[isotope].second;
          DoubleReal intensity = exact_intensities[peptide][isotope];

          if ((intensity / intensity_mono) / (averagineIntensity / averagineIntensity_mono) > model_deviation || (intensity / intensity_mono) / (averagineIntensity / averagineIntensity_mono) < 1 / model_deviation)
          {
            // MISSING PEAK EXCEPTION
            // A missing intensity is allowed if (1) the user allowed it, (2) one of the two peaks is the last isotopic peak of a SILAC peptide and (3) it hasn't occured before.
            if (SILACFiltering::allow_missing_peaks && (isotope == isotopes_per_peptide - 1) && (!missing_peak_seen_yet))
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


    //---------------------------------------------------------------
    // ALL FILTERS PASSED => CREATE DATAPOINT
    //---------------------------------------------------------------
    DataPoint newElement;    // Raw data point at this particular RT and m/z passed all filters. Store it for further clustering.
    newElement.feature_id = SILACFiltering::feature_id;
    newElement.rt = rt;
    newElement.mz = mz;
    newElement.charge = charge;
    newElement.isotopes_per_peptide = isotopes_per_peptide;
    newElement.intensities.insert(newElement.intensities.begin(), exact_intensities.begin(), exact_intensities.end());
    newElement.mass_shifts.insert(newElement.mass_shifts.begin(), mz_peptide_separations.begin(), mz_peptide_separations.end());
    elements.push_back(newElement);

    return true;
  }


  DoubleReal SILACFilter::computeActualMzShift(DoubleReal mz, DoubleReal expectedMzShift, DoubleReal maxMzDeviation)
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
      DoubleReal starting_offset = std::min(mz - SILACFiltering::mz_min, maxMzDeviation);

      // calculate akima interpolation for region around mz (+- maxMzDeviation) and store in vector akimaMz
      // starting position: mz - maxMzDeviation
      // ending position: mz + maxMzDeviation
      // 19 steps, stepwidth = maxMzDeviation / 18
      // akima interpolation for position x
      Size i = 0;
      for (DoubleReal x = mz - starting_offset; x <= mz + maxMzDeviation; x += stepwidth)
      {
        akimaMz[i] = gsl_spline_eval_deriv2 (SILACFiltering::spline_aki, x, SILACFiltering::current_aki);
        ++i;
      }

      // calculate akima interpolation for region around mz + expectedMzShift (- 2 * maxMzDeviation, + maxMzDeviation) and store in vector akimaMzShift
      // starting position: mz - 2 * maxMzDeviation
      // ending position: mz + maxMzDeviation
      // 19 steps, maxMzDeviation / 18 each
      // akima interpolation for position (x + expectedMzShift)
      std::vector<DoubleReal> akimaMzShift(akimaMz_size, 0.0);
      i = 0;
      for (DoubleReal x = mz - starting_offset - maxMzDeviation; x <= mz + maxMzDeviation; x += stepwidth)
      {
        akimaMzShift[i] = gsl_spline_eval_deriv2 (SILACFiltering::spline_aki, x + expectedMzShift, SILACFiltering::current_aki);
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


      for (DoubleReal current_position = 0.0; current_position <= maxMzDeviation; current_position += stepwidth)
      {
        // search for first maximum in + direction and check preconditions
        DoubleReal last_intensity = gsl_spline_eval (spline_correlation, expectedMzShift + current_position - stepwidth, acc_correlation);
        DoubleReal current_intensity = gsl_spline_eval (spline_correlation, expectedMzShift + current_position, acc_correlation);
        DoubleReal next_intensity = gsl_spline_eval (spline_correlation, expectedMzShift + current_position + stepwidth, acc_correlation);

        // search for a current m/z shift larger than the expected one
        // conditions: intensity at position (mz + expectedMzShift + current_position) > intesity_cutoff (intensity calculated with akima interpolation based on "intensities_vec" from SILACFiltering)
        if (current_intensity > last_intensity && current_intensity > next_intensity && gsl_spline_eval (SILACFiltering::spline_aki, mz + expectedMzShift + current_position, SILACFiltering::current_aki) > SILACFiltering::intensity_cutoff) // Why fixed intensity cutoffs?
        {
          gsl_spline_free(spline_correlation);      // free interpolation object
          gsl_interp_accel_free(acc_correlation);     // free accelerator object
          return expectedMzShift + current_position;      // return exact position
        }

        // search for first maximum in - direction and check preconditions
        last_intensity = gsl_spline_eval(spline_correlation, expectedMzShift - current_position - stepwidth, acc_correlation);
        current_intensity = gsl_spline_eval(spline_correlation, expectedMzShift - current_position, acc_correlation);
        next_intensity = gsl_spline_eval(spline_correlation, expectedMzShift - current_position + stepwidth, acc_correlation);

        // search for an current m/z shift smaller than the expected one
        // conditions: intensity at position (mz + expectedMzShift - current_position) > intesity_cutoff (intensity calculated with akima interpolation based on "intensities_vec" from SILACFiltering)
        if (current_intensity > last_intensity && current_intensity > next_intensity && gsl_spline_eval (SILACFiltering::spline_aki, mz + expectedMzShift - current_position, SILACFiltering::current_aki) > SILACFiltering::intensity_cutoff)
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


  DoubleReal SILACFilter::getPeakWidth(DoubleReal mz)
  {
    return 5*(1.889e-7*pow(mz,1.5));
  }

  Int SILACFilter::getSILACType()
  {
    return mz_peptide_separations.size();
  }

  std::vector<DoubleReal> SILACFilter::getPeakPositions()
	{
		peak_positions.clear();
		for (Int peptide = 0; peptide <= numberOfPeptides; peptide++)
		{
			for (Int isotope = 0; isotope < isotopes_per_peptide; isotope++)
			{
				peak_positions.push_back(current_mz + expected_shifts[peptide][isotope]);
			}
		}
		return peak_positions;
  }
	
	std::vector<DoubleReal> SILACFilter::getExpectedMZshifts()
	{
		return expectedMZshifts;
  }

  std::vector<DataPoint> SILACFilter::getElements()
  {
    return elements;
  }

  Int SILACFilter::getCharge()
  {
    return charge;
  }

  std::vector<DoubleReal> SILACFilter::getMassSeparations()
  {
    return mass_separations;
  }

}
