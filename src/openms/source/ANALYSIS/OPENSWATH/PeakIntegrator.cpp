// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

  PeakIntegrator::~PeakIntegrator() = default;

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
    params.setValidStrings("integration_type", {"intensity_sum","simpson","trapezoid"});

    params.setValue("baseline_type", BASELINE_TYPE_BASETOBASE, "The baseline type to use in estimateBackground() based on the peak boundaries. A rectangular baseline shape is computed based either on the minimal intensity of the peak boundaries, the maximum intensity or the average intensity (base_to_base).");
    params.setValidStrings("baseline_type", {"base_to_base","vertical_division","vertical_division_min","vertical_division_max"});

    params.setValue("fit_EMG", "false", "Fit the chromatogram/spectrum to the EMG peak model.");
    params.setValidStrings("fit_EMG", {"false","true"});
  }

  void PeakIntegrator::updateMembers_()
  {
    integration_type_ = (String)param_.getValue("integration_type").toString();
    baseline_type_ = (String)param_.getValue("baseline_type").toString();
    fit_EMG_ = param_.getValue("fit_EMG").toBool();
  }

}
