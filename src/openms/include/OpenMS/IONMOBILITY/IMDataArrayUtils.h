// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eugen Netz, Chris Bielow, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/METADATA/DataArrays.h>

namespace OpenMS
{
  /**
    @brief Interpret ion-mobility array metadata independently of frame conversion.

    Uses the existing PSI-MS vocabulary and vendor-name fallbacks. The shared CV
    files remain a runtime requirement. IMDataConverter retains compatibility
    entry points for both operations.
  */
  class OPENMS_DLLAPI IMDataArrayUtils
  {
  public:
    /**
      @brief Set the array name corresponding to a supported ion-mobility unit.
      @param[in,out] fda Array whose name is updated.
      @param[in] unit MILLISECOND or VSSC.
      @throws Exception::InvalidValue for other units.
    */
    static void setIMUnit(DataArrays::FloatDataArray& fda, const DriftTimeUnit unit);

    /**
      @brief Recognize an ion-mobility array and determine its unit.
      @param[in] fda Array whose name is interpreted.
      @param[out] unit Recognized unit, or NONE for a CV term without usable units.
      Unchanged if the array is not recognized as ion-mobility data.
      @return True if the array describes ion-mobility data.
    */
    static bool getIMUnit(const DataArrays::FloatDataArray& fda, DriftTimeUnit& unit);
  };
}
