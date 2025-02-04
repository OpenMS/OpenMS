// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Juliane Schmachtenberg, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/QC/QCBase.h>
#include <OpenMS/QC/RTAlignment.h>
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
    if (any_of(vdp.begin(), vdp.end(), [](const DataProcessing& dp) {
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
      feature.setMetaValue("rt_align", trafo.apply(feature.getRT()));
      feature.setMetaValue("rt_raw", feature.getRT());
      feature.setMetaValue("rt_align_start", trafo.apply(feature.getConvexHull().getBoundingBox().minX()));
      feature.setMetaValue("rt_align_end", trafo.apply(feature.getConvexHull().getBoundingBox().maxX()));
      feature.setMetaValue("rt_raw_start", feature.getConvexHull().getBoundingBox().minX());
      feature.setMetaValue("rt_raw_end", feature.getConvexHull().getBoundingBox().maxX());
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
  QCBase::Status RTAlignment::requirements() const
  {
    return QCBase::Status() | QCBase::Requires::TRAFOALIGN | QCBase::Requires::POSTFDRFEAT;
  }
} // namespace OpenMS
