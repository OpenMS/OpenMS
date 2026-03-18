// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>

#include <vector>

namespace OpenMS
{
/**
 * @brief fundamental data structure for SplineInterpolatedPeaks
 *
 * In many cases, data points in MS spectra (or chromatograms) are not equidistant in m/z (or RT)
 * but consist of packages of data points separated by wide m/z (or RT) ranges with zero intensity.
 * SplinePackage contains the spline fit of a single set of such data points.
 *
 * @see SplineInterpolatedPeaks
 */
class OPENMS_DLLAPI SplinePackage
{
public:
/**
 * @brief constructor
 */
SplinePackage(std::vector<double> pos, const std::vector<double>& intensity);

/**
 * @brief destructor
 */
~SplinePackage();

/**
 * @brief returns the minimum position for which the spline fit is valid
 */
double getPosMin() const;

/**
 * @brief returns the maximum position for which the spline fit is valid
 */
double getPosMax() const;

/**
 * @brief returns a sensible position step width for the package
 */
double getPosStepWidth() const;

/**
 * @brief returns true if position in [posMin:posMax] interval else false
 */
bool isInPackage(double pos) const;

/**
 * @brief returns interpolated intensity position `pos`
 */
double eval(double pos) const;

private:
/**
 * @brief position limits of the package in the raw data spectrum
 */
double pos_min_;
double pos_max_;

/**
 * @brief sensible position step width with which to scan through the package
 * 
 * @note The step width is rescaled individually in each navigator.
 * @see SplineInterpolatedPeaks::Navigator::getNextPos()
 */
double pos_step_width_;

/**
 * @brief spline object for interpolation of intensity profile
 */
CubicSpline2d spline_;

};

}

