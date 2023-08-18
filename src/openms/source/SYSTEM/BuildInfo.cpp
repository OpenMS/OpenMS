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
//

#include <OpenMS/SYSTEM/BuildInfo.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>


#include <simde/simde-arch.h>

using namespace std;

namespace OpenMS
{

  String Internal::OpenMSOSInfo::getActiveSIMDExtensions()
  {
    StringList ret;
#ifdef SIMDE_ARCH_ARM_NEON
    ret.push_back("neon");
#endif
#ifdef SIMDE_ARCH_X86_SSE
    ret.push_back("SSE");
#endif
#ifdef SIMDE_ARCH_X86_SSE2
    ret.push_back("SSE2");
#endif
#ifdef SIMDE_ARCH_X86_SSE3
    ret.push_back("SSE3");
#endif
#ifdef SIMDE_ARCH_X86_SSE4_1
    ret.push_back("SSE4.1");
#endif
#ifdef SIMDE_ARCH_X86_SSE4_2
    ret.push_back("SSE4.2");
#endif
#ifdef SIMDE_ARCH_X86_AVX
    ret.push_back("AVX");
#endif
#ifdef SIMDE_ARCH_X86_AVX2
    ret.push_back("AVX2");
#endif
#ifdef SIMDE_ARCH_X86_FMA
    ret.push_back("FMA");
#endif

    return ListUtils::concatenate(ret, ", ");
  }

} // namespace OpenMS