// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// This file is NOT part of the upstream Percolator source tree. It lives in
// the vendored percolator/ directory because it adapts upstream code and
// because reviewers syncing from upstream should see it side-by-side with
// SetHandler.cpp (from which it is derived). It is NOT in whitelist.txt and
// must NOT be copied over by sync-from-upstream.sh.

#pragma once

#include <cstddef>
#include <vector>

namespace OpenMS { namespace Internal { namespace Percolator {

/// @brief Reservoir-sample a training subset of row indices.
///
/// Derived from @c SetHandler::readPSMs (upstream @c SetHandler.cpp:234-278
/// at @c rel-3-08-01, commit @c febeef34). The in-memory wrapper over
/// Percolator cannot call @c SetHandler::readPSMs directly because that
/// function reads from a @c std::istream and builds @c PSMDescription
/// objects inline; the in-process path has the data in memory and registers
/// PSMs through @c DataSet::registerPsm instead. This helper replicates
/// only the reservoir-sampling *decision* (which row indices to keep),
/// leaving PSM construction to the caller.
///
/// Algorithm, verbatim from upstream:
/// - Each unique ScanId = (scan_number, exp_mass) gets a single random
///   priority, shared across all rows with the same ScanId. This preserves
///   the "same-spectrum PSMs travel together" invariant.
/// - Random priorities come from @c PseudoRandom::lcg_rand() — the same
///   generator upstream uses; the caller must have seeded it via
///   @c PseudoRandom::setSeed(seed) for reproducibility.
/// - A max-heap of size @c max_train holds (priority, row_idx). A row is
///   pushed iff the heap is not full or its priority is below the current
///   @c upper_limit. On overflow, @c upper_limit is set to the dropped
///   priority (NOT the new top), matching upstream's admission heuristic.
///
/// Returns the sampled indices in ascending row order (deterministic).
/// When @c max_train is 0 or >= @c n_rows, returns @c [0, n_rows) without
/// consulting the RNG.
///
/// @par Re-sync notes
/// When re-vendoring upstream Percolator (see @c sync-from-upstream.sh):
///   1. Open upstream's @c SetHandler.cpp and locate the reservoir block
///      (currently lines 234-278 in @c rel-3-08-01).
///   2. Diff against this function. If the algorithm changed, port the
///      change here.
///   3. The subprocess-parity class tests
///      (@c Percolator_subprocess_parity_test) will flag drift if this
///      function diverges from the shipped @c percolator binary.
///
/// @param scan_numbers length @c n_rows, each row's ScanId.first.
/// @param exp_masses   length @c n_rows, each row's ScanId.second.
/// @param max_train    cap on the returned subset size. 0 means no cap.
/// @return row indices, ascending, of the sampled training subset.
std::vector<std::size_t> sampleTrainingRowIndices(
    const std::vector<int>& scan_numbers,
    const std::vector<double>& exp_masses,
    std::size_t max_train);

}}}  // namespace OpenMS::Internal::Percolator
