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

#ifndef OPENMS_ANALYSIS_OPENSWATH_PEAKINTEGRATOR_H
#define OPENMS_ANALYSIS_OPENSWATH_PEAKINTEGRATOR_H

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

namespace OpenMS
{
  class OPENMS_DLLAPI PeakIntegrator :
    public DefaultParamHandler
  {
public:
    PeakIntegrator();
    virtual ~PeakIntegrator();

    void estimateBackground(
      const MSChromatogram& chromatogram,
      const double& left,
      const double& right
    );

    void integratePeak(
      const MSChromatogram& chromatogram,
      const double& left,
      const double& right
    );

    // internal structure to represent various peak shape metrics
    struct PeakShapeMetrics_ {
      double width_at_5 = 0.0;
      double width_at_10 = 0.0;
      double width_at_50 = 0.0;
      double start_time_at_10 = 0.0;
      double start_time_at_5 = 0.0;
      double start_time_at_50 = 0.0;
      double end_time_at_10 = 0.0;
      double end_time_at_5 = 0.0;
      double end_time_at_50 = 0.0;
      double total_width = 0.0;
      /**
        The tailing factor is a measure of peak tailing.
        It is defined as the distance from the front slope of the peak to the back slope
        divided by twice the distance from the center line of the peak to the front slope,
        with all measurements made at 5% of the maximum peak height.
        tailing_factor = Tf = W0.05/2a
        where W0.05 is peak width at 5% max peak height
        a = min width to peak maximum at 5% max peak height
        b = max width to peak maximum at 5% max peak height
        0.9 < Tf < 1.2
        front Tf < 0.9
        tailing Tf > 1.2
      */
      double tailing_factor = 0.0;
      /**
        The asymmetry factor is a measure of peak tailing.
        It is defined as the distance from the center line of the peak to the back slope
        divided by the distance from the center line of the peak to the front slope,
        with all measurements made at 10% of the maximum peak height.
        asymmetry_factor = As = b/a
        where a is min width to peak maximum at 10% max peak height
        b is max width to peak maximum at 10% max peak height
      */
      double asymmetry_factor = 0.0;
      /**
        The change in baseline divided by the height is
        a way of comparing the influence of the change of baseline on the peak height.
      */
      double baseline_delta_2_height = 0.0;
      /**
        The slope of the baseline is a measure of slope change.
        It is approximated as the difference in baselines between the peak start and peak end.
      */
      double slope_of_baseline = 0.0;
      int points_across_baseline = 0;
      int points_across_half_height = 0;
    };

    void calculatePeakShapeMetrics(
      const MSChromatogram& chromatogram,
      const double& left,
      const double& right,
      PeakShapeMetrics_& peakShapeMetrics
    );

    double getPeakArea() const;
    double getPeakHeight() const;
    double getPeakApexRT() const;
    double getBackgroundHeight() const;
    double getBackgroundArea() const;

    void getDefaultParameters(Param& params);

protected:
    void updateMembers_();

private:
    // parameters
    String integration_type_; // intensity_sum, trapezoid, simpson
    String baseline_type_; // vertical_division, base_to_base
    String peak_model_; // none

    // outputs
    double peak_area_ = 0.0;
    double peak_height_ = -1.0;
    double peak_apex_rt_ = -1.0;
    double background_height_ = 0.0;
    double background_area_ = 0.0;

    // helper
    double simpson(
      MSChromatogram::ConstIterator it_begin,
      MSChromatogram::ConstIterator it_end
    ) const;
  };
}

#endif // OPENMS_ANALYSIS_OPENSWATH_PEAKINTEGRATOR_H
