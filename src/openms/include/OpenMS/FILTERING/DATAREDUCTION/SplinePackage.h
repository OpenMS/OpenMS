// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
#include <algorithm>
#include <iostream>

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
SplinePackage(std::vector<double> pos, std::vector<double> intensity);

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
 * @brief returns interpolated intensity @ position pos
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

