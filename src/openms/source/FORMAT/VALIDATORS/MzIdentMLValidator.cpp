// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/VALIDATORS/MzIdentMLValidator.h>

#include <OpenMS/FORMAT/ControlledVocabulary.h>

using namespace xercesc;
using namespace std;

namespace OpenMS::Internal
{
  MzIdentMLValidator::MzIdentMLValidator(const CVMappings & mapping, const ControlledVocabulary & cv) :
    SemanticValidator(mapping, cv)
  {
    setCheckUnits(true);
  }

  MzIdentMLValidator::~MzIdentMLValidator() = default;
} // namespace OpenMS  // namespace Internal
