// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

namespace OpenMS
{
  /**
    @brief LOWESS (locally weighted scatterplot smoothing).

    A smoothing technique that a quadratic model to localized subsets of the
    data, point by point.

    This is particularly useful for smoothing intensities in spectra or
    chromatograms. In this case, the window size for the smoothing
    should be set proportional to the peak width (see LowessSmoothing parameters).

    Note that this should work best for few datapoints that have strong
    non-linear behavior. For large datasets with mostly linear behavior, use
    FastLowessSmoothing

    @htmlinclude OpenMS_LowessSmoothing.parameters

    @ingroup SignalProcessing
  */
  class OPENMS_DLLAPI LowessSmoothing :
    public DefaultParamHandler
  {
public:
    /// Default constructor
    LowessSmoothing();

    /// Destructor
    ~LowessSmoothing() override;

    typedef std::vector<double> DoubleVector;

    /// Smoothing method that receives x and y coordinates (e.g., RT and intensities) and computes smoothed intensities.
    void smoothData(const DoubleVector &, const DoubleVector &, DoubleVector &);

protected:
    void updateMembers_() override;

private:
    double window_size_;

    double tricube_(double, double);
  };


} // namespace OpenMS
