// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Marc Sturm, Andreas Bertsch, Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmPEPIons.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ConsensusIDAlgorithmPEPIons, "$Id$")

/////////////////////////////////////////////////////////////

ConsensusIDAlgorithm* ptr = nullptr;
ConsensusIDAlgorithm* null_pointer = nullptr;
START_SECTION(ConsensusIDAlgorithmPEPIons())
{
  ptr = new ConsensusIDAlgorithmPEPIons();
  TEST_NOT_EQUAL(ptr, null_pointer);
}
END_SECTION

START_SECTION(~ConsensusIDAlgorithmPEPIons())
{
  delete(ptr);
}
END_SECTION


START_SECTION(void apply(std::vector<PeptideIdentification>& ids))
{
  NOT_TESTABLE // tested by ConsensusID TOPP tool tests
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
