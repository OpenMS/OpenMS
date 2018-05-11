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
// $Author: Erhan Kenar $
// $Maintainer: Timo Sachsenberg$
// --------------------------------------------------------------------------
//

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>

using namespace std;

namespace OpenMS
{
  PeakPickerHiRes::PeakPickerHiRes() :
    DefaultParamHandler("PeakPickerHiRes"),
    ProgressLogger()
  {
    // set default parameter values
    defaults_.setValue("signal_to_noise", 0.0, "Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)");
    defaults_.setMinFloat("signal_to_noise", 0.0);

    defaults_.setValue("spacing_difference_gap", 4.0, "The extension of a peak is stopped if the spacing between two subsequent data points exceeds 'spacing_difference_gap * min_spacing'. 'min_spacing' is the smaller of the two spacings from the peak apex to its two neighboring points. '0' to disable the constraint. Not applicable to chromatograms.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("spacing_difference_gap", 0.0);

    defaults_.setValue("spacing_difference", 1.5, "Maximum allowed difference between points during peak extension, in multiples of the minimal difference between the peak apex and its two neighboring points. If this difference is exceeded a missing point is assumed (see parameter 'missing'). A higher value implies a less stringent peak definition, since individual signals within the peak are allowed to be further apart. '0' to disable the constraint. Not applicable to chromatograms.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("spacing_difference", 0.0);

    defaults_.setValue("missing", 1, "Maximum number of missing points allowed when extending a peak to the left or to the right. A missing data point occurs if the spacing between two subsequent data points exceeds 'spacing_difference * min_spacing'. 'min_spacing' is the smaller of the two spacings from the peak apex to its two neighboring points. Not applicable to chromatograms.", ListUtils::create<String>("advanced"));
    defaults_.setMinInt("missing", 0);

    defaults_.setValue("ms_levels", ListUtils::create<Int>(""), "List of MS levels for which the peak picking is applied. If empty, auto mode is enabled, all peaks which aren't picked yet will get picked. Other scans are copied to the output without changes.");
    defaults_.setMinInt("ms_levels", 1);

    defaults_.setValue("report_FWHM", "false", "Add metadata for FWHM (as floatDataArray named 'FWHM' or 'FWHM_ppm', depending on param 'report_FWHM_unit') for each picked peak.");
    defaults_.setValidStrings("report_FWHM", ListUtils::create<String>("true,false"));
    defaults_.setValue("report_FWHM_unit", "relative", "Unit of FWHM. Either absolute in the unit of input, e.g. 'm/z' for spectra, or relative as ppm (only sensible for spectra, not chromatograms).");
    defaults_.setValidStrings("report_FWHM_unit", ListUtils::create<String>("relative,absolute"));

    // parameters for STN estimator
    defaults_.insert("SignalToNoise:", SignalToNoiseEstimatorMedian< MSSpectrum >().getDefaults());

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }

  PeakPickerHiRes::~PeakPickerHiRes()
  {
  }

  void PeakPickerHiRes::updateMembers_()
  {
    signal_to_noise_ = param_.getValue("signal_to_noise");
    spacing_difference_gap_ = param_.getValue("spacing_difference_gap");
    if (spacing_difference_gap_ == 0.0) spacing_difference_gap_ = std::numeric_limits<double>::infinity();
    spacing_difference_ = param_.getValue("spacing_difference");
    if (spacing_difference_ == 0.0) spacing_difference_ = std::numeric_limits<double>::infinity();
    missing_ = param_.getValue("missing");

    ms_levels_ = getParameters().getValue("ms_levels");
    report_FWHM_ = getParameters().getValue("report_FWHM").toBool();
    report_FWHM_as_ppm_ = getParameters().getValue("report_FWHM_unit")!="absolute";
  }

}
