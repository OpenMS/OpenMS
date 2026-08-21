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
    @brief Reading and writing the charge notation used in fragment ion names

    Producers disagree on whether PeptideHit::PeakAnnotation::annotation spells the charge out:
    TheoreticalSpectrumGenerator writes "y3+", TheoreticalSpectrumGeneratorXLMS writes
    "[alpha|ci$y3]" and leaves the charge to PeakAnnotation::charge, and mzPAF names use "^2".
    A consumer that wants to display a complete ion name therefore has to ask what a given name
    already says, rather than assume. These functions are that one place to ask.

    @ingroup Chemistry
  */
  namespace IonNaming
  {
    /**
      @brief Largest charge magnitude still written as a run of signs

      Beyond this the charge is written as a sign followed by the number, which keeps names readable
      and bounds the allocation for a charge that came from a file.
    */
    constexpr int MAX_REPEATED_SIGNS = 8;

    /**
      @brief The charge notation to append to an ion name, e.g. "+" for 1, "--" for -2

      Up to @p MAX_REPEATED_SIGNS the charge is written as a run of '+' (or '-' in negative mode),
      which is what TheoreticalSpectrumGenerator has always written and what TOPPView has always
      shown. Beyond that the sign is followed by the number ("+12"). That bounds the allocation:
      @p charge can come from a file, where an unbounded run of signs would try to allocate
      gigabytes for a large value.
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

      Understands the three notations that occur in OpenMS: a trailing run of '+'/'-' signs ("y3++"),
      a trailing sign followed by the number ("y3+3"), and the mzPAF caret ("y3-H2O^2", which may be
      followed by mzPAF's mass delta or confidence field, e.g. "y3^2/3.2ppm").
      Only the first line is inspected; further lines of a peak label are free text.

      This is deliberately conservative: a name that does not clearly spell a charge returns 0 rather
      than guessing, so callers can tell "no charge in the name" from "charge 1".
    */
    inline int chargeFromName(const std::string& ion_name)
    {
      const std::string name = ion_name.substr(0, ion_name.find_first_of("\r\n"));
      if (name.empty()) { return 0; }

      // parse name[begin, end) as a charge magnitude; false if it is not one that fits in an int
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
      // mzPAF's following fields ('/' mass delta, '*' confidence). Checked before the trailing-number
      // form below, because mzPAF's mass delta is itself a signed number at the end of the name
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
      @brief The fragment ordinal an ion name spells out ("y12++" -> 12), or 0 if it does not spell one

      The ordinal is the run of digits directly after the ion letter, so it stops at whatever follows --
      a neutral loss, a charge suffix, or an mzPAF field. Deriving it this way rather than by deleting
      the parts that are not the ordinal means a name the caller does not recognise yields 0 instead of
      a wrong number: "y3-H2O1+" has ordinal 3, not 31.

      Returns 0 for names whose ordinal is not in that position, e.g. the cross-link names
      ("[alpha|ci$y3]++"). Callers must range-check the result against the peptide before indexing.
    */
    inline UInt ordinalFromName(const std::string& ion_name)
    {
      // the ordinal follows an ion letter, so anything else in front of the digits is not one:
      // "+3" is a charge and "12" is not an ion name at all
      if (ion_name.empty() || std::isalpha((unsigned char)ion_name[0]) == 0) { return 0; }
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
