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
        @brief Get all unique FAIMS compensation voltages (CVs) that occur in a PeakMap

        - Scans all spectra in the experiment and collects CVs from spectra whose drift time
          unit is DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE.
        - If the data does not contain any FAIMS spectra, an empty set will be returned.
        - The sentinel IMTypes::DRIFTTIME_NOT_SET is ignored; a warning is logged if encountered.

        @param exp The input experiment
        @return Unique FAIMS compensation voltages (in volts)
      */
      static std::set<double> getCompensationVoltages(const PeakMap& exp);
    };

} //end namespace OpenMS
