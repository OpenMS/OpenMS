// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Sven Nahnsen$
// $Authors: Sven Nahnsen$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/IDDecoyProbability.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(IDDecoyProbability, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IDDecoyProbability* ptr = 0;
IDDecoyProbability* null_ptr = 0;
START_SECTION(IDDecoyProbability())
{
	ptr = new IDDecoyProbability();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IDDecoyProbability())
{
	delete ptr;
}
END_SECTION

START_SECTION((IDDecoyProbability(const IDDecoyProbability &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((virtual ~IDDecoyProbability()))
{
  // TODO
}
END_SECTION

START_SECTION((IDDecoyProbability& operator=(const IDDecoyProbability &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((void apply(std::vector< PeptideIdentification > &prob_ids, const std::vector< PeptideIdentification > &fwd_ids, const std::vector< PeptideIdentification > &rev_ids)))
{
  // TODO
}
END_SECTION

START_SECTION((void apply(std::vector< PeptideIdentification > &ids)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



