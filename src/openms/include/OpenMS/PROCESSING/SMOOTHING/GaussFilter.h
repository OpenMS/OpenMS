// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/PROCESSING/SMOOTHING/GaussFilterAlgorithm.h>
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
