// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg, Chris Bielow $
// $Authors: Marc Sturm, Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/StringUtilsSimple.h>

#include <OpenMS/SYSTEM/SIMDe.h>

namespace OpenMS
{

  const char* StringUtils::skipWhitespace(const char* p, const char* p_end)
  {
    const simde__m128i w0 = simde_mm_set1_epi8(' ');
    const simde__m128i w1 = simde_mm_set1_epi8('\t');
    const simde__m128i w2 = simde_mm_set1_epi8('\n');
    const simde__m128i w3 = simde_mm_set1_epi8('\r');

    for (; p <= p_end - 16; p += 16)
    {
      const simde__m128i s = simde_mm_loadu_si128(reinterpret_cast<const simde__m128i*>(p));
      simde__m128i x = simde_mm_cmpeq_epi8(s, w0);
      x = simde_mm_or_si128(x, simde_mm_cmpeq_epi8(s, w1));
      simde__m128i y = simde_mm_cmpeq_epi8(s, w2);
      y = simde_mm_or_si128(y, simde_mm_cmpeq_epi8(s, w3));
      x = simde_mm_or_si128(x, y);
      // invert (i.e any non-spaces will be '1') and convert to a 16-bit int
      // (do not try to convert first and then invert -- not the same!)
      auto non_space = static_cast<uint16_t>(~simde_mm_movemask_epi8(x)); // 16 bit is paramount here. Do not use 32!
      if (non_space != 0)
      {           // some characters are non-whitespace
  #ifdef _MSC_VER // Find the index of first non-whitespace
        unsigned long offset;
        _BitScanForward(&offset, non_space);
        return p + offset;
  #else
        return p + __builtin_ffs(non_space) - 1;
  #endif
      }
    }
    // the remainder
    while (p != p_end)
    {
      if (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t')
        ++p;
      else
        return p;
    }

    return p_end;
  }

  const char* StringUtils::skipNonWhitespace(const char* p, const char* p_end)
  {
    const simde__m128i w0 = simde_mm_set1_epi8(' ');
    const simde__m128i w1 = simde_mm_set1_epi8('\t');
    const simde__m128i w2 = simde_mm_set1_epi8('\n');
    const simde__m128i w3 = simde_mm_set1_epi8('\r');

    for (; p <= p_end - 16; p += 16)
    {
      const simde__m128i s = simde_mm_loadu_si128(reinterpret_cast<const simde__m128i*>(p));
      simde__m128i x = simde_mm_cmpeq_epi8(s, w0);
      x = simde_mm_or_si128(x, simde_mm_cmpeq_epi8(s, w1));
      simde__m128i y = simde_mm_cmpeq_epi8(s, w2);
      y = simde_mm_or_si128(y, simde_mm_cmpeq_epi8(s, w3));
      x = simde_mm_or_si128(x, y);
      // convert to a 16-bit int (i.e any spaces will be '1')
      // (do not try to convert first and then invert -- not the same!)
      auto spaces = static_cast<uint16_t>(simde_mm_movemask_epi8(x)); // 16 bit is paramount here. Do not use 32!
      if (spaces != 0)
      {           // some characters are whitespace
  #ifdef _MSC_VER // Find the index of first whitespace
        unsigned long offset;
        _BitScanForward(&offset, spaces);
        return p + offset;
  #else
        return p + __builtin_ffs(spaces) - 1;
  #endif
      }
    }
    // the remainder
    while (p != p_end)
    {
      if (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t')
        return p;
      else
        ++p;
    }

    return p_end;
  }

} // namespace OpenMS
