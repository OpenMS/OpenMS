// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Athenea Valls Bleher $
// $Authors: Athenea Valls Bleher $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/config.h>
#include <vector>
#include <string>

/**
  @brief Estimates mass error from polysiloxane contaminant peaks

  The algorithm detects polysiloxane contaminant peaks around m/z 445.12003
  and estimates mass error bias and scatter from their ppm deviations.
*/
namespace OpenMS 
{

struct MarkerPeak
{
    double mz;
    double rt;
    double intensity;
};

class OPENMS_DLLAPI MassErrorEstimator
{
    public: 

        /**
        @brief Detect polysiloxane candidate peaks

        Selects candidate polysiloxane peaks based on m/z and intensity filtering.

        @param experiment Input MS experiment
        @return Candidate polysiloxane peaks
        */
        std::vector<MarkerPeak> findPolysiloxaneCandidates (const MSExperiment& experiment); 

        /**
        @brief Select representative peaks across retention time

        Splits candidate peaks into retention time bins and selects
        the most intense peak per bin.

        @param candidates Candidate marker peaks
        @return Selected marker peaks
        */
        std::vector<MarkerPeak> RTdistribution (const std::vector<MarkerPeak>& candidates); // sort candidates and fit distribution and check for unique values

        /**
        @brief Calculate mean ppm mass error

        Computes the mean ppm deviation of selected peaks relative to a reference m/z value.

        @param selected Selected marker peaks
        @param ref Reference m/z value
        @return Mean ppm mass error
        */

        double calculateMean (const std::vector<MarkerPeak>& selected, double ref); 

        /**
        @brief Calculate standard deviation of ppm mass error

        Computes the standard deviation of ppm deviations relative
        to a reference m/z value.

        @param selected Selected marker peaks
        @param ref Reference m/z value
        @return Standard deviation of ppm mass error
        */
        double calculateSD (const std::vector<MarkerPeak>& selected, double ref);

        /**
        @brief Estimate mass error from polysiloxane peaks

        Estimates mass error bias and scatter of an MS experiment
        from detected polysiloxane contaminant peaks.

        @param experiment Input MS experiment
        */
        void estimate (const MSExperiment& experiment);
      

    private: 
    double reference_mz_ = 445.12003; ///< Reference polysiloxane m/z
    double ppm_tolerance_ = 30.0; ///< PPM tolerance window 
    Size intensity_factor_ = 5; ///< Required intensity factor above noise median
    double rt_bins_ = 100.0; ///< Number of retention time bins
    Size min_unique_mz_ = 10; ///< Minimum number of unique m/z values
};
}