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
#include <array>
#include <ostream>

namespace OpenMS
{
  namespace
  {
    using ArrayTerm = IMDataArrayUtils::ArrayTerm;

    // Child terms of 'MS:1002893 ! ion mobility array' in share/OpenMS/CV/psi-ms.obo (data-version 4.2.2).
    // The unit follows the term's has_units relations, checked in this order: MS:1002814 (volt-second per
    // square centimeter) is VSSC, UO:0000028 (millisecond) is MILLISECOND, UO:0000324 (square angstrom) is CCS.
    // The terms with milliseconds also allow seconds (UO:0000010); they keep MILLISECOND.
    // Compiled in, so interpreting the arrays of a spectrum does not depend on the CV files at run time.
    // IMDataArrayUtils_test compares the table with the bundled ontology: update both together.
    constexpr std::array<ArrayTerm, 9> array_terms{{
      {"MS:1002477", "mean ion mobility drift time array", DriftTimeUnit::MILLISECOND},
      {"MS:1002816", "mean ion mobility array", DriftTimeUnit::MILLISECOND},
      {"MS:1003006", "mean inverse reduced ion mobility array", DriftTimeUnit::VSSC},
      {"MS:1003007", "raw ion mobility array", DriftTimeUnit::MILLISECOND},
      {"MS:1003008", "raw inverse reduced ion mobility array", DriftTimeUnit::VSSC},
      {"MS:1003153", "raw ion mobility drift time array", DriftTimeUnit::MILLISECOND},
      {"MS:1003154", "deconvoluted ion mobility array", DriftTimeUnit::MILLISECOND},
      {"MS:1003155", "deconvoluted inverse reduced ion mobility array", DriftTimeUnit::VSSC},
      {"MS:1003156", "deconvoluted ion mobility drift time array", DriftTimeUnit::MILLISECOND},
    }};

    // name of the table term with this accession, empty if there is none
    constexpr std::string_view termName(std::string_view accession)
    {
      for (const ArrayTerm& term : array_terms)
      {
        if (term.accession == accession)
        {
          return term.name;
        }
      }
      return {};
    }

    // array names written by setIMUnit()
    constexpr std::string_view millisecond_array_name = termName("MS:1002816"); // mean ion mobility array
    constexpr std::string_view vssc_array_name = termName("MS:1003008");        // raw inverse reduced ion mobility array
    static_assert(!millisecond_array_name.empty() && !vssc_array_name.empty(), "setIMUnit() uses a term missing from array_terms");
  } // namespace

  std::span<const IMDataArrayUtils::ArrayTerm> IMDataArrayUtils::getArrayTerms()
  {
    return array_terms;
  }

  void IMDataArrayUtils::setIMUnit(DataArrays::FloatDataArray& fda, const DriftTimeUnit unit)
  {
    switch (unit)
    {
      case DriftTimeUnit::MILLISECOND:
        fda.setName(std::string(millisecond_array_name));
        return;
      case DriftTimeUnit::VSSC:
        fda.setName(std::string(vssc_array_name));
        return;
      default:
        // invalid enum ...
        // There is no CV term which can be used to describe the FDA
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unit is not a valid IM unit for float data arrays", driftTimeUnitToString(unit));
    }
  }

  bool IMDataArrayUtils::getIMUnit(const DataArrays::FloatDataArray& fda, DriftTimeUnit& unit)
  {
    // Prefer the PSI-MS names so official names (e.g. MS:1003006
    // "mean inverse reduced ion mobility array") get the CV unit (1/K0 = VSSC),
    // not the UserParam fallback that defaulted them to milliseconds.
    // An exact match against a short table, without allocation or exceptions:
    // this runs for every float array of every spectrum/pixel.
    const std::string& name = fda.getName();
    for (const ArrayTerm& term : array_terms)
    {
      if (term.name == name)
      {
        if (term.unit == DriftTimeUnit::NONE)
        { // a term without units OpenMS can represent (none in the current ontology)
          OPENMS_LOG_WARN << "Warning: FloatDataArray for IonMobility data '" << term.accession << " " << term.name << "' does not contain proper units!" << std::endl;
        }
        unit = term.unit;
        return true;
      }
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
