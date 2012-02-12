// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Steffen Sass, Holger Plattfaut, Bastian Blank $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_DATAREDUCTION_SILACFILTER_H
#define OPENMS_FILTERING_DATAREDUCTION_SILACFILTER_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFiltering.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SILACPattern.h>
#include <OpenMS/FILTERING/DATAREDUCTION/IsotopeDistributionCache.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <queue>
#include <list>

namespace OpenMS
{
  /**
   * @brief Filter to use for SILACFiltering
   *
   * A SILACFilter searches for SILAC patterns, which correspond to the defined mass shifts and charge.
   * Only peaks are taken into account, which were not blacklisted by other filters before e.g. are not part of
   * a SILAC pair yet.
   *
   * @see SILACFiltering
   */   
  class OPENMS_DLLAPI SILACFilter
  {
    private:
      friend class SILACFiltering;

      typedef IsotopeDistributionCache::TheoreticalIsotopePattern TheoreticalIsotopePattern;

      /**
       * @brief mass shift(s) in [Da] to search for
       */
      std::vector<DoubleReal> mass_separations_;

      /**
       * @brief charge of the ions to search for
       */
      Int charge_;

      /**
       * @brief maximal value of which a predicted SILAC feature may deviate from the averagine model
       */
      DoubleReal model_deviation_;

      /**
       * @brief number of peaks per peptide to search for
       */
      Size isotopes_per_peptide_;

      /**
       * @brief minimal intensity of SILAC features
       */
      DoubleReal intensity_cutoff_;

      /**
       * @brief minimal intensity correlation between regions of different peaks
       */
      DoubleReal intensity_correlation_;

      /**
       * @brief flag for missing peaks
       */
      bool allow_missing_peaks_;

      /**
       * Isotope distributions
       */
      static IsotopeDistributionCache *isotope_distribution_;

      /**
       * @brief number of peptides [i.e. number of labelled peptides +1, e.g. for SILAC triplet =3]
       */
      Size number_of_peptides_;

      /**
       * @brief peak positions of SILAC pattern
       */
      std::vector<DoubleReal> peak_positions_;

      /**
       * @brief m/z separtion between individual peptides [e.g. {0 Th, 4 Th, 5 Th}]
       */
      std::vector<DoubleReal> mz_peptide_separations_;

      /**
       * @brief m/z shifts relative to mono-isotopic peak of unlabelled peptide
       */
      std::vector<DoubleReal> expected_mz_shifts_;

      /**
       * @brief distance between isotopic peaks of a peptide in [Th]
       */
      DoubleReal isotope_distance_;

      /**
       * @brief holds the recognized features
       */
      std::vector<SILACPattern> elements_;

      /**
       * @brief m/z at which the filter is currently applied to
       */
      DoubleReal current_mz_;

      /**
       * @brief exact m/z shift of isotopic peaks in a SILAC pattern relative to the mono-isotopic peak of the light peptide, peptides (row) x isotope (column)
       */
      std::vector<std::vector<DoubleReal> > exact_shifts_;

      /**
       * @brief m/z positions mz + exact_shifts in a SILAC pattern, where mz is the m/z of the mono-isotopic peak of light peptide
       */
      std::vector<std::vector<DoubleReal> > exact_mz_positions_;
      
      /**
       * @brief intensities at mz + exact_shifts in a SILAC pattern, where mz is the m/z of the mono-isotopic peak of light peptide
       */
      std::vector<std::vector<DoubleReal> > exact_intensities_;

      /**
       * @brief expected m/z shift of isotopic peaks in a SILAC pattern relative to the mono-isotopic peak of the light peptide, peptides (row) x isotope (column)
       */
      std::vector<std::vector<DoubleReal> > expected_shifts_;

      /**
       * @brief Checks if there exists a SILAC feature at the given position in the raw (interpolated) data, which corresponds to the filter's properties
       * @param rt RT value of the position
       * @param mz m/z value of the position
       */
      bool isSILACPattern_(const MSSpectrum<Peak1D> &, const SILACFiltering::SpectrumInterpolation &, DoubleReal mz, DoubleReal picked_mz, const SILACFiltering &, MSSpectrum<Peak1D> &debug, SILACPattern &pattern);

      /**
       * @brief Checks if there exists a SILAC feature at the given position in the picked data
       */
      bool isSILACPatternPicked_(const MSSpectrum<Peak1D> &, DoubleReal mz, const SILACFiltering &, MSSpectrum<Peak1D> &debug);

      /**
       * @brief Extracts mass shifts and intensities from the raw (interpolated) data
       */
      bool extractMzShiftsAndIntensities_(const MSSpectrum<Peak1D> &, const SILACFiltering::SpectrumInterpolation &, DoubleReal mz, DoubleReal picked_mz, const SILACFiltering &);

      /**
       * @brief Extracts mass shifts and intensities from the picked data
       */
      bool extractMzShiftsAndIntensitiesPicked_(const MSSpectrum<Peak1D> &, DoubleReal mz, const SILACFiltering &);

      /**
       * @brief Extracts mass shifts and intensities from the picked data and returns pattern information
       */
      bool extractMzShiftsAndIntensitiesPickedToPattern_(const MSSpectrum<Peak1D> &, DoubleReal mz, const SILACFiltering &, SILACPattern &pattern);

      /**
       * @brief Checks all peaks against intensity cutoff
       */
      bool intensityFilter_();

      /**
       * @brief Checks peak form correlation between peaks of one isotope
       */
      bool correlationFilter1_(const SILACFiltering::SpectrumInterpolation &, DoubleReal mz, const SILACFiltering &);

      /**
       * @brief Checks peak form correlation between peaks of different isotopes
       */
      bool correlationFilter2_(const SILACFiltering::SpectrumInterpolation &, DoubleReal mz, const SILACFiltering &);

      /**
       * @brief Checks peak intensities against the averagine model
       */
      bool averageneFilter_(DoubleReal mz);

    public:  
      /**
       * @brief detailed constructor for SILAC pair filtering
       * @param mass_separations all mass shifts of the filter
       * @param charge charge of the ions to search for
       * @param model_deviation maximum deviation from the averagine model
       * @param isotopes_per_peptide number of peaks per peptide to search for
       * @param intensity_cutoff ...
       * @param intensity_correlation minimal intensity correlation between regions of different peaks
       * @param allow_missing_peaks flag for missing peaks
       */
      SILACFilter(std::vector<DoubleReal> mass_separations, Int charge, DoubleReal model_deviation, Int isotopes_per_peptide,
          DoubleReal intensity_cutoff, DoubleReal intensity_correlation, bool allow_missing_peaks);

      /**
       * @brief gets the m/z values of all peaks , which belong the last identified feature
       */
      std::vector<DoubleReal> getPeakPositions();

      /**
       * @brief gets the m/z shifts relative to mono-isotopic peak of unlabelled peptide
       */
      const std::vector<DoubleReal>& getExpectedMzShifts();

      /**
       * @brief returns all identified elements
       */
      std::vector<SILACPattern>& getElements();

      /**
       * @brief returns the charge of the filter
       */
      Int getCharge();

      /**
       * @brief returns the mass shifts of the filter in [Da]
       */
      std::vector<DoubleReal>& getMassSeparations();
  };
}

#endif /* SILACFILTER_H_ */
