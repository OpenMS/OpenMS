// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

//
// + Use this header whenever you want to make use of SIMDe, since it ensures
//    + better cross-platform experience (defining missing operators on simde__m128i etc)
//    + you are only using the SIMDe features backed by the build system (i.e. compile flags, e.g. '-mssse3')
// 
// + Include it only in .cpp files, never in a header, since we do not want to expose the
//   SIMDe internals to the outside world
//

// if you want to upgrade to anything more advanced (e.g. SSE4, or AVX),
// add the respective compile flags to <git>/cmake/compiler_flags.cmake
#include <simde/x86/ssse3.h> 

// these operators are defined for GCC/clang, but not in MSVC (TODO: maybe use SFINAE, but that is overkill for the moment)
#ifdef _MSC_VER
inline simde__m128i operator|(const simde__m128i& left, const simde__m128i& right)
{
  return simde_mm_or_si128(left, right);
}
inline simde__m128i& operator|=(simde__m128i& left, const simde__m128i& right)
{
  left = simde_mm_or_si128(left, right);
  return left;
}
inline simde__m128i operator&(const simde__m128i left, const simde__m128i& right)
{
  return simde_mm_and_si128(left, right);
}
#endif

namespace OpenMS
{
}