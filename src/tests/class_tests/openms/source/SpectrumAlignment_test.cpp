// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <iostream>

#include <OpenMS/COMPARISON/SpectrumAlignment.h>
#include <OpenMS/PROCESSING/SCALING/Normalizer.h>
#include <OpenMS/FORMAT/DTAFile.h>

///////////////////////////

START_TEST(SpectrumAlignment, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

SpectrumAlignment* ptr = nullptr;
SpectrumAlignment* nullPointer = nullptr;

START_SECTION(SpectrumAlignment())
	ptr = new SpectrumAlignment();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(virtual ~SpectrumAlignment())
	delete ptr;
END_SECTION

START_SECTION(SpectrumAlignment(const SpectrumAlignment &source))
  SpectrumAlignment sas1;
  Param p(sas1.getParameters());
  p.setValue("tolerance", 0.2);
  sas1.setParameters(p);

  SpectrumAlignment sas2(sas1);

  TEST_EQUAL(sas1.getName(), sas2.getName())
  TEST_EQUAL(sas1.getParameters(), sas2.getParameters())
END_SECTION

START_SECTION(SpectrumAlignment& operator=(const SpectrumAlignment &source))
  SpectrumAlignment sas1;
  Param p(sas1.getParameters());
  p.setValue("tolerance", 0.2);
  sas1.setParameters(p);

  SpectrumAlignment sas2;

	sas2 = sas1;

  TEST_EQUAL(sas1.getName(), sas2.getName())
  TEST_EQUAL(p, sas2.getParameters())

END_SECTION
		    
START_SECTION(template <typename SpectrumType> void getSpectrumAlignment(std::vector< std::pair< Size, Size > > &alignment, const SpectrumType &s1, const SpectrumType &s2) const)
	PeakSpectrum s1, s2;
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s2);

  TOLERANCE_ABSOLUTE(0.01)

	SpectrumAlignment sas1;
	vector<pair<Size, Size > > alignment;

	sas1.getSpectrumAlignment(alignment, s1, s2);
	
	for (vector<pair<Size, Size > >::const_iterator it = alignment.begin(); it != alignment.end(); ++it)
	{
	//	cerr << it->first << " " << it->second << endl;
	}


  TEST_EQUAL(alignment.size(), s1.size())

  s2.resize(100);

	alignment.clear();
  sas1.getSpectrumAlignment(alignment, s1, s2);

  TEST_EQUAL(alignment.size(), 100)

  // TODO: write better test scenario
  PeakSpectrum s3, s4;
  Param p;
  p.setValue("tolerance", 1.01);
  sas1.setParameters(p);

  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("SpectrumAlignment_in1.dta"), s3);
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("SpectrumAlignment_in2.dta"), s4);

  sas1.getSpectrumAlignment(alignment, s3, s4);	

  TEST_EQUAL(alignment.size(), 5)

  ABORT_IF(alignment.size()!=5)
  
  vector<pair<Size, Size > > alignment_result;
  alignment_result.push_back(std::make_pair(0,0));
  alignment_result.push_back(std::make_pair(1,1));
  alignment_result.push_back(std::make_pair(3,3));
  alignment_result.push_back(std::make_pair(4,5));
  alignment_result.push_back(std::make_pair(6,6));

  for (Size i=0;i<5;++i)
  {
    TEST_EQUAL(alignment[i].first, alignment_result[i].first)
    TEST_EQUAL(alignment[i].second, alignment_result[i].second)
  }

  // test relative tolerance alignment
  alignment.clear();
  p.setValue("is_relative_tolerance", "true");
  p.setValue("tolerance", 10.0); // 10 ppm tolerance
  sas1.setParameters(p);
  sas1.getSpectrumAlignment(alignment, s3, s4);
  TEST_EQUAL(alignment.size(), 1)
  ABORT_IF(alignment.size()!=1)
  alignment_result.clear();
  alignment_result.push_back(std::make_pair(6,6));
  for (Size i=0;i<alignment.size();++i)
  {
    TEST_EQUAL(alignment[i].first, alignment_result[i].first)
    TEST_EQUAL(alignment[i].second, alignment_result[i].second)
  }

  alignment.clear();
  p.setValue("is_relative_tolerance", "true");
  p.setValue("tolerance", 1e4); // one percent tolerance
  sas1.setParameters(p);
  sas1.getSpectrumAlignment(alignment, s3, s4);
  for (vector<pair<Size, Size > >::const_iterator it = alignment.begin(); it != alignment.end(); ++it)
  {
    // cerr << it->first << " " << it->second << endl;
  }
  TEST_EQUAL(alignment.size(), 7)
  ABORT_IF(alignment.size()!=7)
  alignment_result.clear();
  alignment_result.push_back(std::make_pair(0,0));
  alignment_result.push_back(std::make_pair(1,1));
  alignment_result.push_back(std::make_pair(2,2));
  alignment_result.push_back(std::make_pair(3,3));
  alignment_result.push_back(std::make_pair(4,5));
  alignment_result.push_back(std::make_pair(5,5));
  alignment_result.push_back(std::make_pair(6,6));
  for (Size i=0;i<7;++i)
  {
    TEST_EQUAL(alignment[i].first, alignment_result[i].first)
    TEST_EQUAL(alignment[i].second, alignment_result[i].second)
  }


END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
