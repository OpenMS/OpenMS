// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#include "PeakTypes.h"

namespace OpenMS {
namespace PipEcho {

  /****************************************************************************/
  bool Acceptor::is_donor_compatible(const Donor& donor, const Window& window) const {
    return donor.feature.getCharge() == this->feature.getCharge() &&
           std::fabs(donor.feature.getRT() - this->feature.getRT()) <= window.rt_tol &&
           std::fabs(donor.feature.getMZ() - this->feature.getMZ()) <= window.mz_tol;
  }
}}
