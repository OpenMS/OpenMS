// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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