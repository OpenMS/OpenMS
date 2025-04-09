// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TARGETED/IncludeExcludeTarget.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(IncludeExcludeTarget, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IncludeExcludeTarget* ptr = nullptr;
IncludeExcludeTarget* null_ptr = nullptr;
START_SECTION(IncludeExcludeTarget())
{
	ptr = new IncludeExcludeTarget();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IncludeExcludeTarget())
{
	delete ptr;
}
END_SECTION

START_SECTION((IncludeExcludeTarget(const IncludeExcludeTarget &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((virtual ~IncludeExcludeTarget()))
{
  // TODO
}
END_SECTION

START_SECTION((void setName(const String &name)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getName() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setPeptideRef(const String &peptide_ref)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getPeptideRef() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setCompoundRef(const String &compound_ref)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getCompoundRef() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setPrecursorMZ(double mz)))
{
  // TODO
}
END_SECTION

START_SECTION((double getPrecursorMZ() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setPrecursorCVTermList(const CVTermList &list)))
{
  // TODO
}
END_SECTION

START_SECTION((void addPrecursorCVTerm(const CVTerm &cv_term)))
{
  // TODO
}
END_SECTION

START_SECTION((const CVTermList& getPrecursorCVTermList() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setProductMZ(double mz)))
{
  // TODO
}
END_SECTION

START_SECTION((double getProductMZ() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setProductCVTermList(const CVTermList &list)))
{
  // TODO
}
END_SECTION

START_SECTION((void addProductCVTerm(const CVTerm &cv_term)))
{
  // TODO
}
END_SECTION

START_SECTION((const CVTermList& getProductCVTermList() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setInterpretations(const std::vector< CVTermList > &interpretations)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<CVTermList>& getInterpretations() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addInterpretation(const CVTermList &interpretation)))
{
  // TODO
}
END_SECTION

START_SECTION((void setConfigurations(const std::vector< Configuration > &configuration)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<Configuration>& getConfigurations() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addConfiguration(const Configuration &configuration)))
{
  // TODO
}
END_SECTION

START_SECTION((void setPrediction(const CVTermList &prediction)))
{
  // TODO
}
END_SECTION

START_SECTION((void addPredictionTerm(const CVTerm &prediction)))
{
  // TODO
}
END_SECTION

START_SECTION((const CVTermList& getPrediction() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setRetentionTime(RetentionTime rt)))
{
  // TODO
}
END_SECTION

START_SECTION((const RetentionTime& getRetentionTime() const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const IncludeExcludeTarget &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator!=(const IncludeExcludeTarget &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((IncludeExcludeTarget& operator=(const IncludeExcludeTarget &rhs)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



