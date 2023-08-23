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