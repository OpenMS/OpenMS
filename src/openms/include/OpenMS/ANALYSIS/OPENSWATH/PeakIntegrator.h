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

    double estimateBackground(
      const MSChromatogram& chromatogram,
      const double& left,
      const double& right
    );

    void integratePeak(
      const MSChromatogram& chromatogram,
      const double& left,
      const double& right
    );

    double getPeakArea() const;
    double getPeakHeight() const;
    double getPeakApexPosition() const;

    void setIntegrationType(const String& integration_type);
    String getIntegrationType() const;

    void setBaselineType(const String& baseline_type);
    String getBaselineType() const;

    void setPeakModel(const String& peak_model);
    String getPeakModel() const;

    void getDefaultParameters(Param& params);

protected:
    void updateMembers_();

private:
    String integration_type_; // intensity_sum, trapezoid, simpson
    String baseline_type_; // vertical_division, base_to_base
    String peak_model_; // none
    double peak_area_ = 0.0;
    double peak_height_ = -1.0;
    double peak_apex_pos_ = -1.0;

    double simpson(
      MSChromatogram::ConstIterator pt_begin,
      MSChromatogram::ConstIterator pt_end
    ) const;
  };
}

#endif // OPENMS_ANALYSIS_OPENSWATH_PEAKINTEGRATOR_H
