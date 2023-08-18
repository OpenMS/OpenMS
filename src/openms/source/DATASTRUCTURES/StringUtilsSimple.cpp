// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2023.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
      y = simde_mm_or_si128(x, simde_mm_cmpeq_epi8(s, w3));
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
    const __m128i w0 = simde_mm_set1_epi8(' ');
    const __m128i w1 = simde_mm_set1_epi8('\t');
    const __m128i w2 = simde_mm_set1_epi8('\n');
    const __m128i w3 = simde_mm_set1_epi8('\r');

    for (; p <= p_end - 16; p += 16)
    {
      const simde__m128i s = simde_mm_loadu_si128(reinterpret_cast<const simde__m128i*>(p));
      simde__m128i x = simde_mm_cmpeq_epi8(s, w0);
      x = simde_mm_or_si128(x, simde_mm_cmpeq_epi8(s, w1));
      simde__m128i y = simde_mm_cmpeq_epi8(s, w2);
      y = simde_mm_or_si128(x, simde_mm_cmpeq_epi8(s, w3));
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
