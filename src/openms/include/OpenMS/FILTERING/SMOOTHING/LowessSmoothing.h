// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------


#ifndef OPENMS_FILTERING_SMOOTHING_LOWESSSMOOTHING_H
#define OPENMS_FILTERING_SMOOTHING_LOWESSSMOOTHING_H

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
#endif // OPENMS_FILTERING_SMOOTHING_LOWESSSMOOTHING_H
