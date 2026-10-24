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

#include <cmath>
#include <cstddef>
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

/////////////////////////////////////////////////////////////
END_TEST
