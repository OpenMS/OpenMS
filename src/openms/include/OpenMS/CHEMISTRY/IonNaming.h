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
#include <climits>
#include <string>

namespace OpenMS
{
  /**
    @brief Charge notation in fragment ion names

    A fragment ion name spells out the charge of the annotated peak in mzPAF's caret notation ("y3^2"),
    see PeptideHit::PeakAnnotation. These functions are the single place where that notation is written
    and read, so producers and consumers cannot drift apart.

    @ingroup Chemistry
  */
  namespace IonNaming
  {
    /**
      @brief The charge notation to append to an ion name, e.g. "^2" for 2, "^-1" for -1

      This is the mzPAF caret notation, which is what MzPAF::toString() writes and what the HUPO-PSI
      peak annotation format specifies. Unlike mzPAF we also spell out charge 1 ("y3^1"): a name
      without a caret means the charge is unknown, which callers have to be able to tell apart from
      "the name says charge 1" (see PeptideHit::PeakAnnotation).
      A charge of 0 (unknown) yields "".
    */
    inline std::string chargeSuffix(int charge)
    {
      if (charge == 0) { return std::string(); }
      return "^" + std::to_string(charge);
    }

    /**
      @brief The charge an ion name spells out, or 0 if it does not spell one

      Reads the caret notation chargeSuffix() writes ("y3-H2O^2", which may be followed by mzPAF's mass
      delta or confidence field, e.g. "y3^2/3.2ppm"), and keeps reading the two notations OpenMS wrote
      before it: a trailing run of '+'/'-' signs ("y3++") and a trailing sign followed by the number
      ("y3+3"). Files carrying those are still out there, so this has to stay.
      Only the first line is inspected; further lines of a peak label are free text.

      This is deliberately conservative: a name that does not clearly spell a charge returns 0 rather than
      guessing, so callers can tell "no charge in the name" from "charge 1".
    */
    inline int chargeFromName(const std::string& ion_name)
    {
      const std::string name = ion_name.substr(0, ion_name.find_first_of("\r\n"));
      if (name.empty()) { return 0; }

      // parse name[begin, end) as a charge magnitude; returns false if it is not one that fits in an int.
      // Widening before the range check matters: the digits come from a file and chargeSuffix() itself can
      // emit 10 of them (INT_MIN -> "-2147483648"), which would overflow a narrower parse.
      auto toCharge = [&name](std::string::size_type begin, std::string::size_type end, bool negative, int& out) {
        if (end <= begin || (end - begin) > 10) { return false; }
        long long value = 0;
        for (std::string::size_type i = begin; i < end; ++i) { value = value * 10 + (name[i] - '0'); }
        if (negative) { value = -value; }
        if (value < (long long)INT_MIN || value > (long long)INT_MAX) { return false; }
        out = (int)value;
        return true;
      };

      // mzPAF: '^' followed by an optionally signed number, then either the end of the name or one of
      // mzPAF's following fields ('/' mass delta, '*' confidence). Checked before the trailing-number form
      // below, because mzPAF's mass delta is itself a signed number at the end of the name
      // ("y4-H2O1^2/-1") and would otherwise be read as the charge.
      const std::string::size_type caret = name.rfind('^');
      if (caret != std::string::npos)
      {
        std::string::size_type c = caret + 1;
        const bool negative = (c < name.size() && name[c] == '-');
        if (c < name.size() && (name[c] == '+' || name[c] == '-')) { ++c; }
        std::string::size_type end = c;
        while (end < name.size() && std::isdigit((unsigned char)name[end])) { ++end; }
        // anything other than a field delimiter after the digits means this is not a charge token
        const bool token_ends = (end == name.size() || name[end] == '/' || name[end] == '*');
        int value = 0;
        if (token_ends && toCharge(c, end, negative, value)) { return value; }
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
        int value = 0;
        if (d > 0 && (name[d - 1] == '+' || name[d - 1] == '-') && toCharge(d, name.size(), name[d - 1] == '-', value))
        {
          return value;
        }
      }
      return 0;
    }

    /**
      @brief The fragment ordinal an ion name spells out ("y12^2" -> 12), or 0 if it does not spell one

      The ordinal is the run of digits directly after the ion letter, so it stops at whatever follows --
      a neutral loss, the charge caret, or an mzPAF field. Deriving it this way rather than by deleting
      the parts that are not the ordinal means a name the caller does not recognise yields 0 instead of a
      wrong number.

      Returns 0 for names whose ordinal is not in that position, e.g. the cross-link names
      ("[alpha|ci$y3]^2") and NuXL's immonium ions ("iY+U-H3PO4^1"). Callers must range-check the result
      against the peptide before indexing with it.
    */
    inline UInt ordinalFromName(const std::string& ion_name)
    {
      std::string::size_type end = 1;
      while (end < ion_name.size() && std::isdigit((unsigned char)ion_name[end])) { ++end; }
      if (end == 1) { return 0; }
      // bounded: an ordinal is at most a peptide length, anything longer is not one
      if (end - 1 > 9) { return 0; }
      UInt value = 0;
      for (std::string::size_type i = 1; i < end; ++i) { value = value * 10 + (UInt)(ion_name[i] - '0'); }
      return value;
    }
  } // namespace IonNaming
} // namespace OpenMS
