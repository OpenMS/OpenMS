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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>

namespace OpenMS
{
  PeakIntegrator::PeakIntegrator() :
    DefaultParamHandler("PeakIntegrator")
  {
    getDefaultParameters(defaults_);
    defaultsToParam_(); // write defaults into Param object param_
  }

  PeakIntegrator::~PeakIntegrator() {}

  PeakIntegrator::PeakArea PeakIntegrator::integratePeak(const MSChromatogram& chromatogram, const double left, const double right) const
  {
    return integratePeak_(chromatogram, left, right);
  }

  PeakIntegrator::PeakArea PeakIntegrator::integratePeak(const MSChromatogram& chromatogram, MSChromatogram::ConstIterator& left, MSChromatogram::ConstIterator& right) const
  {
    return integratePeak_(chromatogram, left->getRT(), right->getRT());
  }

  PeakIntegrator::PeakArea PeakIntegrator::integratePeak(const MSSpectrum& spectrum, const double left, const double right) const
  {
    return integratePeak_(spectrum, left, right);
  }

  PeakIntegrator::PeakArea PeakIntegrator::integratePeak(const MSSpectrum& spectrum, MSSpectrum::ConstIterator& left, MSSpectrum::ConstIterator& right) const
  {
    return integratePeak_(spectrum, left->getMZ(), right->getMZ());
  }

  PeakIntegrator::PeakBackground PeakIntegrator::estimateBackground(const MSChromatogram& chromatogram, const double left, const double right, const double peak_apex_pos) const
  {
    return estimateBackground_(chromatogram, left, right, peak_apex_pos);
  }

  PeakIntegrator::PeakBackground PeakIntegrator::estimateBackground(const MSChromatogram& chromatogram, MSChromatogram::ConstIterator& left, MSChromatogram::ConstIterator& right, const double peak_apex_pos) const
  {
    return estimateBackground_(chromatogram, left->getRT(), right->getRT(), peak_apex_pos);
  }

  PeakIntegrator::PeakBackground PeakIntegrator::estimateBackground(const MSSpectrum& spectrum, const double left, const double right, const double peak_apex_pos) const
  {
    return estimateBackground_(spectrum, left, right, peak_apex_pos);
  }

  PeakIntegrator::PeakBackground PeakIntegrator::estimateBackground(const MSSpectrum& spectrum, MSSpectrum::ConstIterator& left, MSSpectrum::ConstIterator& right, const double peak_apex_pos) const
  {
    return estimateBackground_(spectrum, left->getMZ(), right->getMZ(), peak_apex_pos);
  }

  PeakIntegrator::PeakShapeMetrics PeakIntegrator::calculatePeakShapeMetrics(const MSChromatogram& chromatogram, const double left, const double right, const double peak_height, const double peak_apex_pos) const
  {
    return calculatePeakShapeMetrics_(chromatogram, left, right, peak_height, peak_apex_pos);
  }

  PeakIntegrator::PeakShapeMetrics PeakIntegrator::calculatePeakShapeMetrics(const MSChromatogram& chromatogram, MSChromatogram::ConstIterator& left, MSChromatogram::ConstIterator& right, const double peak_height, const double peak_apex_pos) const
  {
    return calculatePeakShapeMetrics_(chromatogram, left->getRT(), right->getRT(), peak_height, peak_apex_pos);
  }

  PeakIntegrator::PeakShapeMetrics PeakIntegrator::calculatePeakShapeMetrics(const MSSpectrum& spectrum, const double left, const double right, const double peak_height, const double peak_apex_pos) const
  {
    return calculatePeakShapeMetrics_(spectrum, left, right, peak_height, peak_apex_pos);
  }

  PeakIntegrator::PeakShapeMetrics PeakIntegrator::calculatePeakShapeMetrics(const MSSpectrum& spectrum, MSSpectrum::ConstIterator& left, MSSpectrum::ConstIterator& right, const double peak_height, const double peak_apex_pos) const
  {
    return calculatePeakShapeMetrics_(spectrum, left->getMZ(), right->getMZ(), peak_height, peak_apex_pos);
  }

  void PeakIntegrator::getDefaultParameters(Param& params)
  {
    params.clear();

    params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM, "The integration technique to use in integratePeak() and estimateBackground() which uses either the summed intensity, integration by Simpson's rule or trapezoidal integration.");
    params.setValidStrings("integration_type", ListUtils::create<String>("intensity_sum,simpson,trapezoid"));

    params.setValue("baseline_type", BASELINE_TYPE_BASETOBASE, "The baseline type to use in estimateBackground() based on the peak boundaries. A rectangular baseline shape is computed based either on the minimal intensity of the peak boundaries, the maximum intensity or the average intensity (base_to_base).");
    params.setValidStrings("baseline_type", ListUtils::create<String>("base_to_base,vertical_division,vertical_division_min,vertical_division_max"));

    params.setValue("fit_EMG", "false", "Fit the chromatogram/spectrum to the EMG peak model.");
    params.setValidStrings("fit_EMG", ListUtils::create<String>("false,true"));
  }

  void PeakIntegrator::updateMembers_()
  {
    integration_type_ = (String)param_.getValue("integration_type");
    baseline_type_ = (String)param_.getValue("baseline_type");
    fit_EMG_ = param_.getValue("fit_EMG").toBool();
  }

}
