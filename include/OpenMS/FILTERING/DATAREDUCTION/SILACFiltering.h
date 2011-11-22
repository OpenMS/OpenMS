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
// $Authors: Steffen Sass, Holger Plattfaut, Bastian Blank $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_DATAREDUCTION_SILACFILTERING_H
#define OPENMS_FILTERING_DATAREDUCTION_SILACFILTERING_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <list>
#include <map>
#include <vector>

namespace OpenMS
{
  class SILACFilter;

  /**
    @brief Filtering for SILAC data.

    This filtering can be used to extract SILAC features from an MS experiment.
    Several SILACFilters can be added to the filtering to search for specific SILAC patterns.

    @see SILACFilter
  */
  class OPENMS_DLLAPI SILACFiltering
    : public ProgressLogger
  {
    public:
      typedef std::vector<SILACFilter> Filters;

      /**
       * @brief holds all filters used in the filtering
       */
      Filters filters_;

      /**
       * @brief Wrapper class for spectrum interpolation
       */
      class OPENMS_DLLAPI SpectrumInterpolation
      {
        private:
          gsl_interp_accel *current_;
          gsl_spline *spline_;

        public:
          SpectrumInterpolation(const MSSpectrum<> &, const SILACFiltering &);
          ~SpectrumInterpolation();

          DoubleReal operator()(DoubleReal mz) const
          {
            return gsl_spline_eval(spline_, mz, current_);
          }
      };

    private:
      /**
       * @brief minimal intensity of SILAC features
       */
      DoubleReal intensity_cutoff_;

      /**
       * @brief raw data
       */
      MSExperiment<Peak1D>& exp_;

      /**
       * @brief picked data
       */
      MSExperiment<Peak1D> picked_exp_;

      /**
       * @brief picked data seeds
       */
      MSExperiment<Peak1D> picked_exp_seeds_;

      /**
       * Filename base for debugging output
       */
      const String debug_filebase_;

      /**
       * @brief pick data seeds
       */
      void pickSeeds_();

      /**
       * @brief apply filtering to picked data seeds
       */
      void filterSeeds_();

    public:

      /// peak-width equation
      const PeakWidthEstimator::Result peak_width;

      /**
       * @brief detailed constructor
       * @param exp raw data
       * @param intensity_cutoff minimal intensity of SILAC features
       * @param intensity_correlation minimal intensity correlation between regions of different peaks
       * @param allow_missing_peaks flag for missing peaks
       */
      SILACFiltering(MSExperiment<Peak1D>& exp, const PeakWidthEstimator::Result &, const DoubleReal intensity_cutoff, const String debug_filebase_ = "");

      /**
       * @brief adds a new filter to the filtering
       * @param filter filter to add
       */
      void addFilter(SILACFilter& filter);

      /**
       * @brief starts the filtering based on the added filters
       */
      void filterDataPoints();

      /**
       * @brief structure for blacklist
       * @param range m/z and RT interval to be blacklisted
       * @param charge charge of the generating filter
       * @param mass_separations mass separations of the generating filter
       * @param relative_peak_position m/z position of the blacklisted area relative to the mono-isotopic peak of the unlabelled peptide
       */
      struct BlacklistEntry
      {
        DRange<2> range;
        Int charge;
        std::vector<DoubleReal> mass_separations;
        DoubleReal relative_peak_position;
      };

      /**
       * @brief holds the range that is blacklisted for other filters and the filter that generated the blacklist entry
       */
      std::multimap<DoubleReal, BlacklistEntry> blacklist;
  };
}

#endif /* SILACFILTERING_H_ */
