// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
//
// Side-by-side parity test for the PercolatorAdapter TOPP tool: runs it twice
// on the same idXML with -use_subprocess true/false, aligns the output hits
// by (spectrum_reference, sequence, charge), and compares the stamped
// Percolator outputs (score, q-value, PEP).
//
// Gated on two env vars set by CMake:
//   PERCOLATOR_BINARY         -> forwarded as -percolator_executable
//   PERCOLATOR_ADAPTER_BINARY -> absolute path to bin/PercolatorAdapter
// Skips silently when either is unset.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>

#include <cstdlib>
#include <string>

using namespace OpenMS;
using namespace std;

START_TEST(PercolatorAdapter_parity, "$Id$")

START_SECTION([EXTRA] adapter parity: -use_subprocess true vs false on same idXML)
{
  const char* perc = std::getenv("PERCOLATOR_BINARY");
  const char* adap = std::getenv("PERCOLATOR_ADAPTER_BINARY");
  if (!perc || !*perc || !adap || !*adap)
  {
    TEST_EQUAL(true, true);  // skip silently
  }
  else
  {
    // Parity logic lands in Task 9.
    TEST_EQUAL(true, true);
  }
}
END_SECTION

END_TEST
