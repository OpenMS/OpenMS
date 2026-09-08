// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/IONMOBILITY/IMDataArrayUtils.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <ostream>

namespace OpenMS
{
  void IMDataArrayUtils::setIMUnit(DataArrays::FloatDataArray& fda, const DriftTimeUnit unit)
  {
    const auto& cv = ControlledVocabulary::getPSIMSCV();
    switch (unit)
    {
      case DriftTimeUnit::MILLISECOND: 
        fda.setName(cv.getTerm("MS:1002816").name); // MS:1002816 ! mean ion mobility array
        return;
      case DriftTimeUnit::VSSC:
        fda.setName(cv.getTerm("MS:1003008").name); // MS:1003008 ! raw inverse reduced ion mobility array
        return;
      default:
        // invalid enum ...
        // There is no CV term which can be used to describe the FDA
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unit is not a valid IM unit for float data arrays", driftTimeUnitToString(unit));
    }
  }

  bool IMDataArrayUtils::getIMUnit(const DataArrays::FloatDataArray& fda, DriftTimeUnit& unit)
  {
    const auto& cv = ControlledVocabulary::getPSIMSCV();

    // Prefer PSI-MS ontology lookup so official names (e.g. MS:1003006
    // "mean inverse reduced ion mobility array") get the CV unit (1/K0 = VSSC),
    // not the UserParam fallback that defaulted them to milliseconds.
    // Unknown names yield nullptr, so non-CV arrays cost no exception here —
    // this runs for every float array of every spectrum/pixel.
    const ControlledVocabulary::CVTerm* cv_term = cv.checkAndGetTermByName(fda.getName());
    if (cv_term != nullptr && cv.isChildOf(cv_term->id, "MS:1002893")) // child of generic 'ion mobility array'?
    {
      if (cv_term->units.contains("MS:1002814"))
      { // MS:1002814 ! volt-second per square centimeter
        unit = DriftTimeUnit::VSSC;
      }
      else if (cv_term->units.contains("UO:0000028"))
      { // UO:0000028 ! millisecond
        unit = DriftTimeUnit::MILLISECOND;
      }
      else if (cv_term->units.contains("UO:0000324"))
      { // UO:0000324 ! square angstrom (CCS)
        unit = DriftTimeUnit::CCS;
      }
      else
      { // fallback
        OPENMS_LOG_WARN << "Warning: FloatDataArray for IonMobility data '" << cv_term->id << " " << cv_term->name << "' does not contain proper units!" << std::endl;
        unit = DriftTimeUnit::NONE;
      }
      return true;
    }

    // Fallbacks for non-standard / vendor UserParam names
    // (Mobi-DIK "Ion Mobility", MSConvert "inverse reduced ion mobility", …).
    if (StringUtils::hasPrefix(fda.getName(), Constants::UserParam::MEAN_INVERSE_REDUCED_ION_MOBILITY_ARRAY) ||
        StringUtils::hasPrefix(fda.getName(), Constants::UserParam::INVERSE_REDUCED_ION_MOBILITY))
    {
      // These names denote 1/K0 (Vs/cm^2), not drift time in ms.
      unit = DriftTimeUnit::VSSC;
      return true;
    }
    if (StringUtils::hasPrefix(fda.getName(), Constants::UserParam::ION_MOBILITY))
    {
      if (StringUtils::hasSubstring(fda.getName(), "MS:1002815") || StringUtils::hasSubstring(fda.getName(), "MS:1003006"))
      {
        unit = DriftTimeUnit::VSSC;
      }
      else if (StringUtils::hasSubstring(fda.getName(), "MS:1002954"))
      {
        unit = DriftTimeUnit::CCS;
      }
      else
      {
        unit = DriftTimeUnit::MILLISECOND;
      }
      return true;
    }
    return false;
  }

} // namespace OpenMS
