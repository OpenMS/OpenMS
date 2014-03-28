// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Erhan Kenar$
// --------------------------------------------------------------------------
//

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>

#include <vector>

using namespace std;

namespace OpenMS
{
  PeakPickerHiRes::PeakPickerHiRes() :
    DefaultParamHandler("PeakPickerHiRes"),
    ProgressLogger()
  {
    // set default parameter values
    defaults_.setValue("signal_to_noise", 1.0, "Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)");
    defaults_.setMinFloat("signal_to_noise", 0.0);

    defaults_.setValue("spacing_difference", 1.5, "Maximum allowed distance between peaks in multiples of the minimal difference. A higher value is implies a less stringent peak definition since individual signals within the peaks are allowed to further apart. E.g. if the value is set to 1.5 and in a peak the minimal spacing between peaks is 10 mDa, then only signals at most 15 mDa apart will be added to the peak.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("spacing_difference", std::numeric_limits<unsigned>::min()); //must be > 0

    defaults_.setValue("ms1_only", "false", "If true, peak picking is only applied to MS1 scans. Other scans are copied to the output without changes.");
    defaults_.setValidStrings("ms1_only", ListUtils::create<String>("true,false"));

    // parameters for SNT estimator
    defaults_.insert("SignalToNoise:", SignalToNoiseEstimatorMedian< MSSpectrum<Peak1D> >().getDefaults());

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
    spacing_difference_ = param_.getValue("spacing_difference");
  }

}
