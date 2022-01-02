// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilterAlgorithm.h>
#include <OpenMS/KERNEL/StandardTypes.h>

namespace OpenMS
{
  /**
    @brief This class represents a Gaussian lowpass-filter which works on uniform as well as on non-uniform profile data.

    Gaussian filters are important in many signal processing,
    image processing, and communication applications. These filters are characterized by narrow bandwidths,
    sharp cutoffs, and low passband ripple. A key feature of Gaussian filters is that the Fourier transform of a
    Gaussian is also a Gaussian, so the filter has the same response shape in both the time and frequency domains.
    The coefficients \f$ \emph{coeffs} \f$ of the Gaussian-window with length \f$ \emph{frameSize} \f$ are calculated
    from the gaussian distribution
    \f[ \emph{coeff}(x) = \frac{1}{\sigma \sqrt{2\pi}} e^{\frac{-x^2}{2\sigma^2}} \f]
    where \f$ x=[-\frac{frameSize}{2},...,\frac{frameSize}{2}] \f$ represents the window area and \f$ \sigma \f$
    is the standard derivation.

    @note The wider the kernel width the smoother the signal (the more detail information get lost!).
          Use a Gaussian filter kernel which has approximately the same width as your mass peaks,
          whereas the Gaussian peak width corresponds approximately to 8*sigma.

        @note The data must be sorted according to ascending m/z!

        @htmlinclude OpenMS_GaussFilter.parameters

    @ingroup SignalProcessing
  */
//#define DEBUG_FILTERING

  class OPENMS_DLLAPI GaussFilter :
    public ProgressLogger,
    public DefaultParamHandler
  {
public:
    /// Constructor
    GaussFilter();

    /// Destructor
    ~GaussFilter() override = default;

    /**
      @brief Smoothes an MSSpectrum containing profile data.

      Convolutes the filter and the profile data and writes the result back to the spectrum.

      @exception Exception::IllegalArgument is thrown, if the @em gaussian_width parameter is too small.
    */
    void filter(MSSpectrum & spectrum);

    void filter(MSChromatogram & chromatogram);

    /**
      @brief Smoothes an MSExperiment containing profile data.

      @exception Exception::IllegalArgument is thrown, if the @em gaussian_width parameter is too small.
    */
    void filterExperiment(PeakMap & map);

protected:

    GaussFilterAlgorithm gauss_algo_;

    /// The spacing of the pre-tabulated kernel coefficients
    double spacing_;

    bool write_log_messages_ = false;

    // Docu in base class
    void updateMembers_() override;
  };

} // namespace OpenMS
