// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/QuantitativeExperimentalDesign.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(QuantitativeExperimentalDesign, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

QuantitativeExperimentalDesign* ptr = nullptr;
QuantitativeExperimentalDesign* null_ptr = nullptr;
START_SECTION(QuantitativeExperimentalDesign())
{
	ptr = new QuantitativeExperimentalDesign();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~QuantitativeExperimentalDesign())
{
	delete ptr;
}
END_SECTION

START_SECTION((virtual ~QuantitativeExperimentalDesign()))
{
  // TODO
}
END_SECTION

START_SECTION((void applyDesign2Resolver(ProteinResolver &resolver, TextFile &file, StringList &fileNames)))
{
  // TODO
}
END_SECTION

START_SECTION((void applyDesign2Quantifier(PeptideAndProteinQuant &quantifier, TextFile &file, StringList &fileNames)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



