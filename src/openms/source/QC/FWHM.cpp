// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/QC/FWHM.h>

namespace OpenMS
{
  void FWHM::compute(FeatureMap& features)
  {
    for (auto& f : features)
    {
      if (f.metaValueExists("FWHM")) // from FF-Centroided
      {
        for (auto& pi : f.getPeptideIdentifications())
        {
          pi.setMetaValue("FWHM", f.getMetaValue("FWHM"));
        }
      }
      else if (f.metaValueExists("model_FWHM")) // from FF-Identification
      {
        for (auto& pi : f.getPeptideIdentifications())
        {
          pi.setMetaValue("FWHM", f.getMetaValue("model_FWHM")); // use 'FWHM' as target to make the name unique for downstream processing
        }
      }
      else
      {
        // throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Metavalue 'FWHM' or 'model_FWHM' is missing for a feature in a FeatureMap. Please check your FeatureFinder
        // reports FWHM using these metavalues or add a new mapping here.");
      }
    }
  }

  const String& FWHM::getName() const
  {
    static const String& name = "FWHM";
    return name;
  }

  QCBase::Status FWHM::requirements() const
  {
    return QCBase::Status() | QCBase::Requires::POSTFDRFEAT;
  }
} // namespace OpenMS
