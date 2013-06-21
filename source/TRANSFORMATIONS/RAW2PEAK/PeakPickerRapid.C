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
// $Maintainer: Erhan Kenar $
// $Authors: Erhan Kenar $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerRapid.h>

#include <vector>

using namespace std;

namespace OpenMS
{
PeakPickerRapid::PeakPickerRapid()
    : DefaultParamHandler("PeakPickerRapid"),
      ProgressLogger()
{
    // set default parameter values
    // defaults_.setValue("signal_to_noise", 1.0, "Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)");
    // defaults_.setMinFloat("signal_to_noise", 0.0);

    defaults_.setValue("intensity_type", "maxpeak", "A peak intensity is stored as peak area or as the maximum peak's intensity.");
    defaults_.setValidStrings("intensity_type", StringList::create("peakarea,maxpeak"));

    defaults_.setValue("ms1_only", "true", "If true, peak picking is only applied to MS1 scans. Other scans are copied to the output without changes.");
    defaults_.setValidStrings("ms1_only", StringList::create("true,false"));

    // parameters for SNTestimator config ...

    // write defaults into Param object param_
    defaultsToParam_();


    // initialize class members
    // signal_to_noise_ = param_.getValue("signal_to_noise");
}

PeakPickerRapid::~PeakPickerRapid()
{
}

void PeakPickerRapid::updateMembers_()
{
    // signal_to_noise_ = param_.getValue("signal_to_noise");
}

} // namespace OpenMS
