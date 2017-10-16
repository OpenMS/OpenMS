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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_SPECTRUMEXTRACTOR_H
#define OPENMS_ANALYSIS_OPENSWATH_SPECTRUMEXTRACTOR_H

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/KERNEL/MSSpectrum.h> // MSSpectrum
#include <OpenMS/DATASTRUCTURES/String.h> // String
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h> // DefaultParamHandler
#include <OpenMS/CONCEPT/LogStream.h> // LOG_DEBUG
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h> // SavitzkyGolayFilter
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h> // GaussFilter
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h> // PeakPickerHiRes

namespace OpenMS
{
  /**
  @brief TODO BRIEF EXPLANATION HERE
  
  TODO MULTILINE
  EXPLANATION
  HERE

  */
  class OPENMS_DLLAPI SpectrumExtractor :
    public DefaultParamHandler
  {
public:
    SpectrumExtractor();
    virtual ~SpectrumExtractor();

    void setRTWindow(const double& rt_window);
    double getRTWindow() const;

    void setMinScore(const double& min_score);
    double getMinScore() const;

    void setMinForwardMatch(const double& min_forward_match);
    double getMinForwardMatch() const;

    void setMinReverseMatch(const double& min_reverse_match);
    double getMinReverseMatch() const;

    void setMZTolerance(const double& mz_tolerance);
    double getMZTolerance() const;

    void setMZToleranceUnits(const String& mz_tolerance_units);
    String getMZToleranceUnits() const;

    void setSGolayFrameLength(const UInt& sgolay_frame_length);
    UInt getSGolayFrameLength() const;

    void setSGolayPolynomialOrder(const UInt& sgolay_polynomial_order);
    UInt getSGolayPolynomialOrder() const;

    void setGaussWidth(const double& gauss_width);
    double getGaussWidth() const;

    void setUseGauss(const bool& use_gauss);
    bool getUseGauss() const;

    void setSignalToNoise(const double& signal_to_noise);
    double getSignalToNoise() const;

    void getDefaultParameters(Param& params);

    void pickSpectrum(const MSSpectrum& spectrum, MSSpectrum& picked_spectrum);

protected:
    /// overridden function from DefaultParamHandler to keep members up to date, when a parameter is changed
    void updateMembers_();

private:
    double rt_window_;
    double min_score_;
    double min_forward_match_;
    double min_reverse_match_;
    double mz_tolerance_;
    String mz_tolerance_units_;

    UInt sgolay_frame_length_;
    UInt sgolay_polynomial_order_;
    double gauss_width_;
    bool use_gauss_;
    double signal_to_noise_;
  };
}

#endif // OPENMS_ANALYSIS_OPENSWATH_SPECTRUMEXTRACTOR_H

