// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/QC/PeptideMass.h>

namespace OpenMS
{
  void PeptideMass::compute(FeatureMap& features)
  {
    features.applyFunctionOnPeptideIDs(
      [](PeptideIdentification& pi) {
        if (pi.getHits().empty())
        {
          return;
        }
        auto& hit = pi.getHits()[0];
        hit.setMetaValue("mass", (pi.getMZ() - Constants::PROTON_MASS_U) * hit.getCharge());
      },
      true);
  }

  const String& PeptideMass::getName() const
  {
    static const String& name = "PeptideMass";
    return name;
  }

  QCBase::Status PeptideMass::requirements() const
  {
    return QCBase::Status() | QCBase::Requires::POSTFDRFEAT;
  }
} // namespace OpenMS
