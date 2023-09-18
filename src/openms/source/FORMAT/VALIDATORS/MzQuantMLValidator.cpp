// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/VALIDATORS/MzQuantMLValidator.h>

using namespace xercesc;
using namespace std;

namespace OpenMS::Internal
{

  MzQuantMLValidator::MzQuantMLValidator(const CVMappings & mapping, const ControlledVocabulary & cv) :
    SemanticValidator(mapping, cv)
  {
    setCheckUnits(true);
  }

  MzQuantMLValidator::~MzQuantMLValidator() = default;

} // namespace OpenMS   // namespace Internal
