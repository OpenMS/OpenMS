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
#include <OpenMS/MATH/MISC/BSpline2d.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>

namespace OpenMS
{
  /**
   * @brief Rough estimation of the peak width at m/z
   *
   * Based on the peaks of the dataset (peak position & width) and the peak
   * boundaries as reported by the PeakPickerHiRes, the typical peak width is
   * estimated for arbitrary m/z using a spline interpolation.
   *
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

