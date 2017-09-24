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
// $Maintainer: Douglas McCloskey $
// $Authors: Douglas McCloskeyt $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>

#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/QcMLFile.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{

  MRMFeatureQC::MRMFeatureQC() :
    DefaultParamHandler("MRMFeatureQC")
  {
    //TODO
    defaults_.setValue("stop_after_feature", -1, "Stop finding after feature (ordered by intensity; -1 means do not stop).");
    defaults_.setValue("stop_after_intensity_ratio", 0.0001, "Stop after reaching intensity ratio");
    defaults_.setValue("min_peak_width", -1.0, "Minimal peak width (s), discard all peaks below this value (-1 means no action).", ListUtils::create<String>("advanced"));

    defaults_.setValue("background_subtraction", "none", "Try to apply a background subtraction to the peak (experimental). The background is estimated at the peak boundaries, either the smoothed or the raw chromatogram data can be used for that.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("background_subtraction", ListUtils::create<String>("none,smoothed,original"));

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }

  MRMFeatureQC::~MRMFeatureQC()
  {
  }

  MRMFeatureQC& MRMFeatureQC::operator=(const MRMFeatureQC& rhs)
  {
    if (&rhs == this)
      return *this;

    // don't copy parameters

    return *this;
  }

  void MRMFeatureQC::updateMembers_()
  {
    //TODO
    stop_after_feature_ = (int)param_.getValue("stop_after_feature");
  }

  void MRMFeatureQC::FilterFeatureMap(FeatureMap& features)
  {
    //TODO
  }
  
  void MRMFeatureQC::FeatureMapToAttachment(FeatureMap& features, QcMLFile::Attachment& attachment)
  {
    //TODO
  }

}

