// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/TVControllerBase.h>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
  TVControllerBase::TVControllerBase(TOPPViewBase* parent):
    tv_(parent)
  {  
  }

  void TVControllerBase::activateBehavior()
  {
    // no special handling of activation is default
  }

  void TVControllerBase::deactivateBehavior()
  {
    // no special handling of deactivation is default
  }
}
