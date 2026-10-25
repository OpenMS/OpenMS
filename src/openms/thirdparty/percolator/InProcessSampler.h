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
#include <cstdint>
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

/// @brief The spectra of an out-of-fold capped training (OpenMS; not upstream).
struct OutOfFoldSample
{
  /// per row: its cross-validation fold (0 to num_folds - 1), the same for all rows of a spectrum
  std::vector<unsigned char> fold;
  /// the rows of the sampled spectra, ascending
  std::vector<std::size_t> sampled_rows;
  /// spectra in the input and among the sampled rows
  std::size_t num_spectra = 0;
  std::size_t num_sampled_spectra = 0;
  /// the rows of the validation spectra (@p validation_rows > 0), ascending; never sampled
  std::vector<std::size_t> validation_rows;
};

/// @brief Out-of-fold capped training: assigns every spectrum of the input to a cross-validation fold and draws a
///        random sample of whole spectra of at most @p max_train rows from the entire input.
///
/// A spectrum is a (spec_file, scan) pair, Percolator's own key of its cross-validation folds (CrossValidation:
/// createXvalSetsBySpectrum), so all PSMs of a spectrum (also with different precursor masses) share a fold and are
/// sampled together. Each spectrum gets a 64-bit priority and a fold from a hash of (@p seed, spec_file, scan): the
/// draw depends neither on the order of the rows nor on the threads. The sample is the spectra of the lowest priorities
/// (ties by spec_file, scan) for as long as their rows fit into @p max_train (it stops at the first spectrum that does
/// not fit): a random-priority prefix of the spectra within a row budget. The inclusion probability depends on a
/// spectrum's size only through the boundary (a spectrum larger than the remaining budget ends the sample, so the
/// sample may hold fewer rows than @p max_train, by less than one spectrum's rows), which matters only when spectra are
/// not small against the cap. The cap is global (one
/// sample for the cross-validation, as Percolator's --subset-max-train): each fold's model trains on the sampled
/// spectra of the other folds (about 2/3 of the sample).
///
/// @param spec_files length n_rows: each row's spectrum file number
/// @param scans      length n_rows: each row's scan number
/// @param max_train  cap on the sampled rows (> 0)
/// @param seed       the seed of the draw (Percolator's seed)
/// With @p validation_rows > 0 (the learning-curve mode), the spectra at the END of the priority order are set aside
/// first, for as long as their rows fit into @p validation_rows (from the last spectrum backwards, stopping at the first
/// that does not fit), and the sample is drawn from the other spectra. The validation spectra do not depend on
/// @p max_train, and the samples of increasing @p max_train are nested (prefixes of the same order).
///
/// @param num_folds  number of cross-validation folds (> 0)
/// @param validation_rows  rows set aside as validation spectra (0 = none)
OutOfFoldSample drawOutOfFoldSample(const std::vector<int>& spec_files,
                                   const std::vector<int>& scans,
                                   std::size_t max_train,
                                   std::uint64_t seed,
                                   unsigned int num_folds,
                                   std::size_t validation_rows = 0);

}}}  // namespace OpenMS::Internal::Percolator
