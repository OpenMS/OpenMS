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

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h>

#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/QcMLFile.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{

  MRMFeatureFilter::MRMFeatureFilter() :
    DefaultParamHandler("MRMFeatureFilter")
  {
    defaults_.setValue("flag_or_filter", "flag", "Flag or Filter (i.e., remove) Components or transitions that do not pass the QC.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("flag_or_filter", ListUtils::create<String>("flag,filter"));
    
    defaults_.setValue("report_xic", "false", "Embed an image of the XIC in the QC report.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("report_xic", ListUtils::create<String>("true,false"));

    defaults_.setValue("report_tic", "false", "Embed an image of the TIC in the QC report.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("report_tic", ListUtils::create<String>("true,false"));

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }

  MRMFeatureFilter::~MRMFeatureFilter()
  {
  }

  MRMFeatureFilter& MRMFeatureFilter::operator=(const MRMFeatureFilter& rhs)
  {
    if (&rhs == this)
      return *this;

    // don't copy parameters

    return *this;
  }

  void MRMFeatureFilter::updateMembers_()
  {
    flag_or_filter_ = (String)param_.getValue("flag_or_filter");
    report_xic_ = (bool)param_.getValue("report_xic");
    report_tic_ = (bool)param_.getValue("report_tic");
  }

  void MRMFeatureFilter::FilterFeatureMap(FeatureMap& features)
  { 
    // initialize the new feature map
    if (flag_or_filter_ == "filter")
    {
      FeatureMap features_filtered;
    }
    
    // initialize QC variables
    std::map<String,MRMFeatureFilterFile>::iterator feature_qc_it;

    // initialize variables
    String component_name; //i.e., transition_id
    String IS_component_name; //i.e., internal standard transition_id
    String component_group_name; //i.e., peptideRef
    double calculated_concentration;
    bool qc_pass;
    String concentration_units;// iterate through each component_group/feature     

    for (size_t feature_it = 0; feature_it < features.size(); ++feature_it)
    {
      component_group_name = (String)features[feature_it].getMetaValue("PeptideRef");

      // iterate through each component/sub-feature
      for (size_t sub_it = 0; sub_it < features[feature_it].getSubordinates().size(); ++sub_it)
      {
        component_name = (String)features[feature_it].getSubordinates()[sub_it].getMetaValue("native_id"); 
        qc_pass = false;

        // iterate through multi-feature/multi-sub-feature QCs/filters
        if (sub_it == 0)
        {
          //TODO
        }


        // iterate through feature/sub-feature QCs/filters

      }
    }
  }
  
  void MRMFeatureFilter::FeatureMapToAttachment(FeatureMap& features, QcMLFile::Attachment& attachment)
  {
    //TODO
  }

}

