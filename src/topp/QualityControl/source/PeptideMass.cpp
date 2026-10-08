// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/FeatureMap.h>
#include "PeptideMass.h"

namespace OpenMS
{
  void PeptideMass::compute(FeatureMap& features)
  {
    // the identifications of the features and the unassigned ones
    QCBase::annotateIdentifications(features, [](Feature*, std::vector<QCBase::AnnotatedIdentification>& identifications) {
      for (auto& id : identifications)
      {
        if (id.top == nullptr)
        {
          continue;
        }
        id.top->setMetaValue("mass", (id.mz - Constants::PROTON_MASS_U) * id.top->getCharge());
      }
    });
  }

  const std::string& PeptideMass::getName() const
  {
    static const std::string& name = "PeptideMass";
    return name;
  }

  QCBase::Status PeptideMass::requirements() const
  {
    return QCBase::Status() | QCBase::Requires::POSTFDRFEAT;
  }
} // namespace OpenMS
