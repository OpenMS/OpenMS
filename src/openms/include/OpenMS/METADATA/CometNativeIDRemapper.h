// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <string>
#include <unordered_map>

namespace OpenMS
{
/**
  @brief Rewrite spectrum native IDs to a monotonic "index=N" form and translate PSM references back.

  Motivation (CometAdapter): Comet's bundled mzParser greps for "scan=" anywhere in a spectrum native
  ID and uses @c atoi of the following digits as a sort key. Bruker DDA native IDs of the form
  "frame=F scan=S precursor=P" carry a non-monotonic "scan=S" (the TIMS isolation-window start), so
  mzParser detects "unsorted" scans and falls into a @c std::sort with a strict-weak-ordering-violating
  comparator -> undefined behavior / segfault. The workaround is to rewrite native IDs to "index=N"
  (which mzParser indexes by counter, not by @c atoi("scan=")) before handing the mzML to Comet, then
  translate the PSM spectrum references on Comet's output back to the original native IDs.

  This namespace exposes that rewrite/restore pair as a small, SDK-free, testable unit so the round-trip
  can be verified without a Comet binary or vendor data.

  @ingroup Metadata
*/
namespace CometNativeIDRemapper
{
  /// MetaValue key under which the original native ID is stashed by rewriteToIndex().
  inline const std::string ORIGINAL_NATIVE_ID = "original_native_id";

  /**
    @brief Rewrite every spectrum's native ID to "index=N" (monotonic, 0-based).

    The original native ID of each spectrum is preserved under the @ref ORIGINAL_NATIVE_ID MetaValue
    so it can be restored later by translateReferencesBack().

    @param exp Experiment whose spectra are rewritten in place.
    @return the number of spectra rewritten.
  */
  inline Size rewriteToIndex(MSExperiment& exp)
  {
    Size idx = 0;
    for (auto& spec : exp.getSpectra())
    {
      spec.setMetaValue(ORIGINAL_NATIVE_ID, spec.getNativeID());
      spec.setNativeID("index=" + StringUtils::toStr(idx++));
    }
    return idx;
  }

  /**
    @brief Translate PSM spectrum references from the rewritten "index=N" form back to the original IDs.

    Restores the spectrum reference of each PeptideIdentification that points at a rewritten "index=N"
    native ID to the original native ID recorded by rewriteToIndex(). References that do not match any
    rewritten spectrum are left untouched. No-op when @p exp is empty or contains no spectrum with an
    @ref ORIGINAL_NATIVE_ID MetaValue.

    @param exp Experiment previously passed to rewriteToIndex() (read-only here).
    @param pids PeptideIdentifications whose spectrum references are translated in place.
  */
  inline void translateReferencesBack(const MSExperiment& exp, PeptideIdentificationList& pids)
  {
    if (exp.empty()) { return; }

    // Build rewritten-id -> original-id map once (O(N) vs O(N*M) linear scan).
    std::unordered_map<std::string, std::string> id_map;
    id_map.reserve(exp.size());
    for (const auto& spec : exp.getSpectra())
    {
      if (spec.metaValueExists(ORIGINAL_NATIVE_ID)) { id_map.emplace(spec.getNativeID(), spec.getMetaValue(ORIGINAL_NATIVE_ID).toString()); }
    }
    if (id_map.empty()) { return; }

    for (auto& pid : pids)
    {
      auto it = id_map.find(pid.getSpectrumReference());
      if (it != id_map.end()) { pid.setSpectrumReference(it->second); }
    }
  }
} // namespace CometNativeIDRemapper
} // namespace OpenMS
