// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>

#include <cctype>
#include <string>

namespace OpenMS
{
  /**
    @brief Charge notation in fragment ion names

    A fragment ion name spells out the charge of the annotated peak, see PeptideHit::PeakAnnotation.
    These two functions are the single place where that notation is written and read, so producers and
    consumers cannot drift apart.

    @ingroup Chemistry
  */
  namespace IonNaming
  {
    /// beyond this magnitude the charge is written as a number instead of a run of signs
    constexpr int MAX_REPEATED_SIGNS = 8;

    /**
      @brief The charge notation to append to an ion name, e.g. "+" for 1, "--" for -2

      Up to @p MAX_REPEATED_SIGNS the charge is written as a run of '+' (or '-' in negative mode), which is
      what TheoreticalSpectrumGenerator has always done. Beyond that the sign is followed by the number
      ("+12"). That keeps names readable and, more importantly, bounds the allocation: @p charge can come
      from a file (the mzIdentML IonType charge attribute, Sage's fragment_charge column), where an
      unbounded run of signs would try to allocate gigabytes for a large value.
      A charge of 0 (unknown) yields "".
    */
    inline std::string chargeSuffix(int charge)
    {
      if (charge == 0) { return std::string(); }
      const char sign = (charge < 0) ? '-' : '+';
      // negate in a wider type: std::abs(INT_MIN) is undefined behaviour
      const long long magnitude = (charge < 0) ? -static_cast<long long>(charge) : static_cast<long long>(charge);
      if (magnitude <= MAX_REPEATED_SIGNS) { return std::string((Size)magnitude, sign); }
      return sign + std::to_string(magnitude);
    }

    /**
      @brief The charge an ion name spells out, or 0 if it does not spell one

      Understands the three notations that occur in OpenMS: a trailing run of '+'/'-' signs ("y3++"), a
      trailing sign followed by the number ("y3+3"), and the mzPAF caret ("y3-H2O^2", which may be followed
      by mzPAF's mass delta or confidence field, e.g. "y3^2/3.2ppm").
      Only the first line is inspected; further lines of a peak label are free text.

      This is deliberately conservative: a name that does not clearly spell a charge returns 0 rather than
      guessing, so callers can tell "no charge in the name" from "charge 1".
    */
    inline int chargeFromName(const std::string& ion_name)
    {
      const std::string name = ion_name.substr(0, ion_name.find_first_of("\r\n"));
      if (name.empty()) { return 0; }

      // a charge never has more digits than this; anything longer is not a charge (and would overflow)
      const std::string::size_type MAX_DIGITS = 9;

      // mzPAF: '^' followed by an optionally signed number, optionally followed by further mzPAF fields.
      // Checked before the trailing-number form below, because mzPAF's mass delta is itself a signed
      // number at the end of the name ("y4-H2O1^2/-1") and would otherwise be read as the charge.
      const std::string::size_type caret = name.rfind('^');
      if (caret != std::string::npos)
      {
        std::string::size_type c = caret + 1;
        const bool negative = (c < name.size() && name[c] == '-');
        if (c < name.size() && (name[c] == '+' || name[c] == '-')) { ++c; }
        std::string::size_type end = c;
        while (end < name.size() && std::isdigit((unsigned char)name[end])) { ++end; }
        if (end > c && (end - c) <= MAX_DIGITS)
        {
          const int value = std::stoi(name.substr(c, end - c));
          return negative ? -value : value;
        }
      }

      const char last = name.back();

      // trailing run of signs, e.g. "y3+", "b2++", "c1-", "a3-B-"
      if (last == '+' || last == '-')
      {
        std::string::size_type run = name.size();
        while (run > 0 && name[run - 1] == last) { --run; }
        const int count = (int)(name.size() - run);
        return (last == '-') ? -count : count;
      }

      // trailing sign plus number, e.g. "y3+3" (note "y1-H2O1" must NOT match: 'O' is not a sign)
      if (std::isdigit((unsigned char)last))
      {
        std::string::size_type d = name.size();
        while (d > 0 && std::isdigit((unsigned char)name[d - 1])) { --d; }
        if (d > 0 && (name[d - 1] == '+' || name[d - 1] == '-') && (name.size() - d) <= MAX_DIGITS)
        {
          const int value = std::stoi(name.substr(d));
          return (name[d - 1] == '-') ? -value : value;
        }
      }
      return 0;
    }
  } // namespace IonNaming
} // namespace OpenMS
