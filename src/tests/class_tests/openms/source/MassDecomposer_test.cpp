// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/MassDecomposer.h>
///////////////////////////

using namespace OpenMS;
using namespace ims;
using namespace std;

START_TEST(MassDecomposer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((virtual ~MassDecomposer()))
{
  // MassDecomposer is an abstract base class, without any implementation
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual bool exist(value_type mass)=0))
{
  // MassDecomposer is an abstract base class, without any implementation
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual decomposition_type getDecomposition(value_type mass)=0))
{
  // MassDecomposer is an abstract base class, without any implementation
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual decompositions_type getAllDecompositions(value_type mass)=0))
{
  // MassDecomposer is an abstract base class, without any implementation
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual decomposition_value_type getNumberOfDecompositions(value_type mass)=0))
{
  // MassDecomposer is an abstract base class, without any implementation
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



