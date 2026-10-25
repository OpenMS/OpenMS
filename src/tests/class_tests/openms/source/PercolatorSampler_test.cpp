// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

// Vendored Percolator (OpenMS::Internal::Percolator): this test links the
// static OpenMS_Percolator library directly, its symbols are hidden in libOpenMS.
#include "InProcessSampler.h"
#include "PseudoRandom.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

using namespace OpenMS::Internal::Percolator;

START_TEST(PercolatorSampler, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION((std::vector<std::size_t> sampleTrainingRowIndices(const std::vector<int>& scan_numbers, const std::vector<double>& exp_masses, std::size_t max_train)))
{
  // The training cap (Percolator subset_max_train, ProSE rescoring:subset_max_train) samples at random from ALL rows,
  // seeded, by spectrum (all PSMs of a spectrum together): not the first max_train rows (user rule, 24 Oct 2026).
  const std::size_t n_spectra = 50000, per_spectrum = 2, n_rows = n_spectra * per_spectrum, max_train = 10000;
  std::vector<int> scans(n_rows);
  std::vector<double> masses(n_rows);
  for (std::size_t r = 0; r < n_rows; ++r)
  {
    scans[r] = static_cast<int>(r / per_spectrum);
    masses[r] = 400.0 + 0.001 * static_cast<double>(r / per_spectrum);
  }
  PseudoRandom::setSeed(1);
  const std::vector<std::size_t> sample = sampleTrainingRowIndices(scans, masses, max_train);
  TEST_TRUE(sample.size() >= max_train - per_spectrum && sample.size() <= max_train + per_spectrum)
  // ascending, and not the prefix [0, max_train)
  bool ascending = true;
  for (std::size_t i = 1; i < sample.size(); ++i) ascending = ascending && sample[i - 1] < sample[i];
  TEST_TRUE(ascending)
  TEST_TRUE(sample.back() >= n_rows - n_rows / 100)
  // spread over the whole input: every tenth of the rows holds about a tenth of the sample
  std::vector<std::size_t> per_decile(10, 0);
  for (std::size_t r : sample) ++per_decile[r * 10 / n_rows];
  const double expected = static_cast<double>(sample.size()) / 10.0;
  for (std::size_t d : per_decile) TEST_TRUE(std::fabs(static_cast<double>(d) - expected) < 0.2 * expected)
  // the PSMs of a spectrum stay together
  std::vector<int> taken(n_spectra, 0);
  for (std::size_t r : sample) ++taken[r / per_spectrum];
  bool together = true;
  for (int t : taken) together = together && (t == 0 || t == static_cast<int>(per_spectrum));
  TEST_TRUE(together)
  // seeded: the same seed draws the same rows, another seed others
  PseudoRandom::setSeed(1);
  TEST_TRUE(sampleTrainingRowIndices(scans, masses, max_train) == sample)
  PseudoRandom::setSeed(2);
  TEST_TRUE(sampleTrainingRowIndices(scans, masses, max_train) != sample)
  // no cap, or a cap of at least the rows: all rows
  TEST_EQUAL(sampleTrainingRowIndices(scans, masses, 0).size(), n_rows)
  TEST_EQUAL(sampleTrainingRowIndices(scans, masses, n_rows).size(), n_rows)
}
END_SECTION

START_SECTION((OutOfFoldSample drawOutOfFoldSample(const std::vector<int>& spec_files, const std::vector<int>& scans, std::size_t max_train, std::uint64_t seed, unsigned int num_folds)))
{
  // Out-of-fold capped training: every spectrum (spec_file, scan) gets a fold; the sample is whole spectra drawn at
  // random from the entire input (not a prefix), seeded, independent of the order of the rows.
  // 2 files x 15000 scans; the rows of a spectrum are not adjacent (three hits: rows s, s + n, s + 2n), with 1 to 3
  // hits per spectrum
  const std::size_t n_spectra = 30000;
  std::vector<int> files, scans;
  for (std::size_t hit = 0; hit < 3; ++hit)
  {
    for (std::size_t s = 0; s < n_spectra; ++s)
    {
      if (hit >= 1 + s % 3) continue;
      files.push_back(static_cast<int>(s % 2));
      scans.push_back(static_cast<int>(s / 2)); // the same scan number in both files: different spectra
    }
  }
  const std::size_t n_rows = scans.size(), max_train = 12000;
  const OutOfFoldSample oof = drawOutOfFoldSample(files, scans, max_train, 1, 3);
  TEST_EQUAL(oof.fold.size(), n_rows)
  TEST_EQUAL(oof.num_spectra, n_spectra)
  TEST_TRUE(oof.sampled_rows.size() <= max_train && oof.sampled_rows.size() + 3 > max_train)
  bool ascending = true;
  for (std::size_t i = 1; i < oof.sampled_rows.size(); ++i) ascending = ascending && oof.sampled_rows[i - 1] < oof.sampled_rows[i];
  TEST_TRUE(ascending)
  // whole spectra: all rows of a spectrum share the fold and the sampling
  std::map<std::pair<int, int>, std::pair<int, int>> fold_and_sampled; // spectrum -> (fold, sampled)
  std::vector<char> sampled(n_rows, 0);
  for (std::size_t r : oof.sampled_rows) sampled[r] = 1;
  bool together = true, in_range = true;
  std::vector<std::size_t> per_fold(3, 0), sampled_per_fold(3, 0);
  for (std::size_t r = 0; r < n_rows; ++r)
  {
    in_range = in_range && oof.fold[r] < 3;
    const std::pair<int, int> v(oof.fold[r], sampled[r]);
    auto it = fold_and_sampled.emplace(std::make_pair(files[r], scans[r]), v).first;
    together = together && it->second == v;
    ++per_fold[oof.fold[r]];
    if (sampled[r]) ++sampled_per_fold[oof.fold[r]];
  }
  TEST_TRUE(together)
  TEST_TRUE(in_range)
  TEST_EQUAL(fold_and_sampled.size(), n_spectra)
  TEST_EQUAL(oof.num_sampled_spectra, static_cast<std::size_t>(std::count_if(fold_and_sampled.begin(), fold_and_sampled.end(),
    [](const auto& e) { return e.second.second == 1; })))
  // the folds and the sample's folds about a third each
  for (std::size_t f = 0; f < 3; ++f)
  {
    TEST_TRUE(std::fabs(static_cast<double>(per_fold[f]) - n_rows / 3.0) < 0.05 * n_rows / 3.0)
    TEST_TRUE(std::fabs(static_cast<double>(sampled_per_fold[f]) - oof.sampled_rows.size() / 3.0) < 0.1 * oof.sampled_rows.size() / 3.0)
  }
  // spread over the whole input: every tenth of each hit block holds about its share of the sample; the last rows too
  std::vector<std::size_t> per_decile(10, 0);
  for (std::size_t r : oof.sampled_rows) ++per_decile[r * 10 / n_rows];
  for (std::size_t d = 0; d < 10; ++d)
  {
    std::size_t rows_in_decile = 0;
    for (std::size_t r = d * n_rows / 10; r < (d + 1) * n_rows / 10; ++r) ++rows_in_decile;
    const double expected = static_cast<double>(rows_in_decile) * oof.sampled_rows.size() / n_rows;
    TEST_TRUE(std::fabs(static_cast<double>(per_decile[d]) - expected) < 0.15 * expected)
  }
  TEST_TRUE(oof.sampled_rows.back() >= n_rows - n_rows / 100)
  // spectra with 1, 2 and 3 PSMs are drawn alike (no preference for small spectra): each third of the spectra about
  // a third of the sampled spectra
  std::vector<std::size_t> by_size(3, 0);
  for (const auto& e : fold_and_sampled)
  {
    if (e.second.second) ++by_size[static_cast<std::size_t>(2 * e.first.second + e.first.first) % 3];
  }
  for (std::size_t k = 0; k < 3; ++k)
  {
    TEST_TRUE(std::fabs(static_cast<double>(by_size[k]) - oof.num_sampled_spectra / 3.0) < 0.1 * oof.num_sampled_spectra / 3.0)
  }
  // seeded: the same seed draws the same, another seed another sample and other folds
  const OutOfFoldSample again = drawOutOfFoldSample(files, scans, max_train, 1, 3);
  TEST_TRUE(again.fold == oof.fold && again.sampled_rows == oof.sampled_rows)
  const OutOfFoldSample other = drawOutOfFoldSample(files, scans, max_train, 2, 3);
  TEST_TRUE(other.fold != oof.fold && other.sampled_rows != oof.sampled_rows)
  // independent of the order of the rows: reversed rows give each spectrum the same fold and sampling
  std::vector<int> rfiles(files.rbegin(), files.rend()), rscans(scans.rbegin(), scans.rend());
  const OutOfFoldSample rev = drawOutOfFoldSample(rfiles, rscans, max_train, 1, 3);
  std::vector<char> rsampled(n_rows, 0);
  for (std::size_t r : rev.sampled_rows) rsampled[r] = 1;
  bool order_free = rev.sampled_rows.size() == oof.sampled_rows.size();
  for (std::size_t r = 0; r < n_rows; ++r)
  {
    order_free = order_free && rev.fold[n_rows - 1 - r] == oof.fold[r] && rsampled[n_rows - 1 - r] == sampled[r];
  }
  TEST_TRUE(order_free)
  // a spectrum larger than the cap is never split: one spectrum of all rows, cap below -> nothing sampled
  const OutOfFoldSample one = drawOutOfFoldSample(std::vector<int>(10, 0), std::vector<int>(10, 5), 9, 1, 3);
  TEST_EQUAL(one.sampled_rows.size(), 0)
  TEST_EQUAL(one.num_spectra, 1)
  TEST_EXCEPTION(std::invalid_argument, drawOutOfFoldSample(files, scans, 0, 1, 3))
  // validation spectra (learning curve): set aside first, the same for every cap, never sampled, whole spectra; the
  // samples of growing caps are nested and leave the folds unchanged
  {
    const OutOfFoldSample small = drawOutOfFoldSample(files, scans, 5000, 1, 3, 8000);
    const OutOfFoldSample large = drawOutOfFoldSample(files, scans, 20000, 1, 3, 8000);
    TEST_TRUE(small.validation_rows == large.validation_rows)
    TEST_TRUE(small.validation_rows.size() <= 8000 && small.validation_rows.size() + 3 > 8000)
    TEST_TRUE(small.fold == oof.fold && large.fold == oof.fold)
    std::vector<char> in_validation(n_rows, 0), in_large(n_rows, 0);
    for (std::size_t r : large.validation_rows) in_validation[r] = 1;
    for (std::size_t r : large.sampled_rows) in_large[r] = 1;
    bool disjoint = true, nested = true, whole = true;
    for (std::size_t r : large.sampled_rows) disjoint = disjoint && !in_validation[r];
    for (std::size_t r : small.sampled_rows) nested = nested && in_large[r];
    std::map<std::pair<int, int>, int> v_of_spectrum;
    for (std::size_t r = 0; r < n_rows; ++r)
    {
      auto it = v_of_spectrum.emplace(std::make_pair(files[r], scans[r]), in_validation[r]).first;
      whole = whole && it->second == in_validation[r];
    }
    TEST_TRUE(disjoint)
    TEST_TRUE(nested)
    TEST_TRUE(whole)
    TEST_TRUE(large.sampled_rows.size() <= 20000 && large.sampled_rows.size() + 3 > 20000)
    // without validation spectra, the sample of a cap is the same as before
    TEST_TRUE(drawOutOfFoldSample(files, scans, max_train, 1, 3, 0).sampled_rows == oof.sampled_rows)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
