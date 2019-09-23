//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Chris Bielow $
// $Authors: Juliane Schmachtenberg, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/QC/RTAlignment.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/QC/QCBase.h>

#include <algorithm>

using namespace std;

namespace OpenMS
{
  // take the original retention time before map alignment and use the transformation information of the post alignment trafoXML
  // for calculation of the post map alignment retention times.
  void RTAlignment::compute(FeatureMap& features, const TransformationDescription& trafo) const
  {
    if (features.empty())
    {
      OPENMS_LOG_WARN << "The FeatureMap is empty.\n";
    }

    // if featureMap after map alignment was handed, return Exception
    auto vdp = features.getDataProcessing(); // get a copy to avoid calling .begin() and .end() on two different temporaries
    if (any_of(vdp.begin(), vdp.end(), [](const DataProcessing& dp){
          return (find(dp.getProcessingActions().begin(), dp.getProcessingActions().end(), DataProcessing::ProcessingAction::ALIGNMENT) != dp.getProcessingActions().end());
        }))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Metric RTAlignment received a featureXML AFTER map alignment, but needs a featureXML BEFORE map alignment!");
    }

    // set meta values for original retention time and aligned retention time (after map alignment)
    for (Feature& feature : features)
    {
      for (PeptideIdentification& peptide_ID : feature.getPeptideIdentifications())
      {
          peptide_ID.setMetaValue("rt_align", trafo.apply(peptide_ID.getRT()));
          peptide_ID.setMetaValue("rt_raw", peptide_ID.getRT());
      }
    }

    // same for unassigned PepIDs
    compute(features.getUnassignedPeptideIdentifications(), trafo);
  }

  void RTAlignment::compute(std::vector<PeptideIdentification>& ids, const TransformationDescription& trafo) const
  {
    // set meta values for all unasssigned PeptideIdentifications
    for (PeptideIdentification& id : ids)
    {
      id.setMetaValue("rt_align", trafo.apply(id.getRT()));
      id.setMetaValue("rt_raw", id.getRT());
    }
  }
  
  const String& RTAlignment::getName() const
  {
    return name_;
  }

  // required input files
  QCBase::Status RTAlignment::requires() const
  {
    return QCBase::Status() | QCBase::Requires::TRAFOALIGN | QCBase::Requires::POSTFDRFEAT;
  }
}
