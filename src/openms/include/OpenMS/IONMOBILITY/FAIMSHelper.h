// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <set>

namespace OpenMS
{

    /**
      @brief Helper functions for FAIMS data

      FAIMSHelper contains convenience functions to deal with FAIMS
      compensation voltages and related data.

    */
    class OPENMS_DLLAPI FAIMSHelper
    {
    public:
      virtual ~FAIMSHelper() {}

      /**
        @brief Get all FAIMS compensation voltages that occur in a PeakMap

        If the data is not FAIMS, an empty set will be returned.

        @param exp The PeakMap with FAIMS data
      */
      static std::set<double> getCompensationVoltages(const PeakMap& exp);
    };

} //end namespace OpenMS
