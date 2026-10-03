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

#include <span>
#include <string_view>

namespace OpenMS
{
  /**
    @brief Interpret ion-mobility array metadata independently of frame conversion.

    Array names are matched against the PSI-MS ion-mobility array terms (getArrayTerms()),
    then against vendor and UserParam names. The terms are compiled in from the psi-ms.obo
    shipped with OpenMS, so no CV file is read at run time. A psi-ms.obo replaced at run
    time does not change which names are recognized: a new ontology term needs an update
    of the table, which IMDataArrayUtils_test compares with the bundled ontology.
    IMDataConverter retains compatibility entry points for both operations.
  */
  class OPENMS_DLLAPI IMDataArrayUtils
  {
  public:
    /// A PSI-MS ion-mobility array term and the unit getIMUnit() assigns to it
    struct ArrayTerm
    {
      std::string_view accession; ///< CV accession, e.g. "MS:1003008"
      std::string_view name;      ///< canonical CV name, as used for the array name
      DriftTimeUnit unit;         ///< unit derived from the term's units in the ontology
    };

    /**
      @brief The PSI-MS ion-mobility array terms that getIMUnit() recognizes.

      These are all child terms of 'MS:1002893 ! ion mobility array' in the bundled psi-ms.obo.
      The unit of a term is the first of VSSC, MILLISECOND and CCS that the term allows, or NONE.
      Terms allowing seconds and milliseconds get MILLISECOND.
    */
    static std::span<const ArrayTerm> getArrayTerms();

    /**
      @brief Set the array name corresponding to a supported ion-mobility unit.

      MILLISECOND names the array 'mean ion mobility array' (MS:1002816), VSSC
      'raw inverse reduced ion mobility array' (MS:1003008).

      @param[in,out] fda Array whose name is updated.
      @param[in] unit MILLISECOND or VSSC.
      @throws Exception::InvalidValue for other units.
    */
    static void setIMUnit(DataArrays::FloatDataArray& fda, const DriftTimeUnit unit);

    /**
      @brief Recognize an ion-mobility array and determine its unit.

      A name equal to the name of a term in getArrayTerms() (case-sensitive) gets that term's unit.
      Otherwise, names starting with 'mean inverse reduced ion mobility array' or 'inverse reduced
      ion mobility' are VSSC. Names starting with 'Ion Mobility' are VSSC if they contain MS:1002815
      or MS:1003006, CCS if they contain MS:1002954, and MILLISECOND otherwise.

      @param[in] fda Array whose name is interpreted.
      @param[out] unit Recognized unit, or NONE for a CV term without usable units.
      Unchanged if the array is not recognized as ion-mobility data.
      @return True if the array describes ion-mobility data.
    */
    static bool getIMUnit(const DataArrays::FloatDataArray& fda, DriftTimeUnit& unit);
  };
}
