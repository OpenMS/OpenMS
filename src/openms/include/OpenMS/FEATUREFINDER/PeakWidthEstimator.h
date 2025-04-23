// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/MATH/MISC/BSpline2d.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>

namespace OpenMS
{
  /**
    @brief Rough estimation of the peak width at m/z
   
    Based on the peaks of the dataset (peak position & width) and the peak
    boundaries as reported by the PeakPickerHiRes, the typical peak width is
    estimated for arbitrary m/z using a spline interpolation.
   
   */
  class OPENMS_DLLAPI PeakWidthEstimator
  {

    public:

    /**
    * @brief constructor
    * 
    * @param exp_picked m/z positions of picked peaks
    * @param boundaries corresponding peak widths
    *
    * @throw Exception::UnableToFit if the B-spline initialisation fails.
    */
    PeakWidthEstimator(const PeakMap & exp_picked, const std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > & boundaries);

    /**
    * @brief returns the estimated peak width at m/z
    * 
    * @note If the mz value is outside the interpolation range, the peak width
    *       value of the maximal/minimal value is used.
    *
    * @throw Exception::InvalidValue if the peak width estimation returns a negative value.
    */
    double getPeakWidth(double mz);

    /**
    * @brief destructor
    */
    virtual ~PeakWidthEstimator();

    private:        

    /// hide default constructor
    PeakWidthEstimator();

    /**
     * @brief B-spline for peak width interpolation
     */
    BSpline2d* bspline_;

    /**
    * @brief m/z range of peak width interpolation
    */
    double mz_min_;
    double mz_max_;
            
  };
}

