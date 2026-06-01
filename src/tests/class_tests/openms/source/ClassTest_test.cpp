// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/PrecisionWrapper.h>

#include <fstream>

using namespace OpenMS;

bool intentionally_failed_tests_okay = true;

/*  This macro turns a previous failure into success.  This is used for
testing the test macros.  It should follow the subtest immediately; preferably
on the same line of code.  */
#define FAILURE_IS_SUCCESS																							\
	if ( !TEST::this_test )																								\
	{																																			\
		TEST::this_test = true;																							\
		TEST::test = intentionally_failed_tests_okay;												\
		if (TEST::verbose > 1)																							\
		{																																		\
			TEST::initialNewline();																						\
			stdcout << __FILE__ ":" <<  __LINE__ <<													\
				": note:  The preceeding test was supposed to fail intentionally.  =>  SUCCESS" << \
				std::endl;																											\
		}																																		\
	}																																			\
	else																																	\
	{																																			\
		TEST::this_test = false;																						\
		intentionally_failed_tests_okay = false;														\
		TEST::test = intentionally_failed_tests_okay;												\
		if (TEST::verbose > 1)																							\
		{																																		\
			TEST::initialNewline();																						\
			stdcout << __FILE__ ":" <<  __LINE__ <<													\
				" error:  The preceeding test was supposed to fail, but it did not.  =>  FAILURE" << \
				std::endl;																											\
		}																																		\
	}

// the only error here is that we do not follow the coding convention ;-)
void
throw_a_Precondition_Exception()
{
  throw OpenMS::Exception::Precondition(__FILE__,__LINE__, OPENMS_PRETTY_FUNCTION,
    "intentional Exception::Preconditon raised by throw_a_Precondition_Exception()");
}

// the only error here is that we do not follow the coding convention ;-)
void
throw_a_Postcondition_Exception()
{
  throw OpenMS::Exception::Postcondition(__FILE__,__LINE__, OPENMS_PRETTY_FUNCTION,
    "intentional Exception::Postconditon raised by throw_a_Postcondition_Exception()");
}


START_TEST(ClassTest, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION("empty section without NOT_TESTABLE")
{
	STATUS("This test should complain about no subtests being performed.");
}
END_SECTION

START_SECTION("empty section with NOT_TESTABLE")
{
	STATUS("This test should NOT complain about no subtests being performed.");
	NOT_TESTABLE;
}
END_SECTION

START_SECTION("TOLERANCE_ABSOLUTE()")
{
	TOLERANCE_ABSOLUTE(0.55);
	TEST_EQUAL(Internal::ClassTest::absdiff_max_allowed, 0.55);
}
END_SECTION

START_SECTION("TOLERANCE_RELATIVE()")
{
	TOLERANCE_RELATIVE(0.66);
	TEST_EQUAL(Internal::ClassTest::ratio_max_allowed, 0.66);
}
END_SECTION

START_SECTION("NEW_TMP_FILE()")
	std::string tmp_filename;
	NEW_TMP_FILE(tmp_filename);
	TEST::this_test = (!tmp_filename.empty());
	TEST_EQUAL(!tmp_filename.empty(), true);
END_SECTION

START_SECTION("TEST_EQUAL()")
{
  TEST_EQUAL(2, 3); FAILURE_IS_SUCCESS;
  TEST_EQUAL(2, 2);
  TEST_EQUAL(2, 3.0); FAILURE_IS_SUCCESS;
  TEST_EQUAL(2, 2.0);
  TEST_EQUAL(2, 2.0f);
  TEST_EQUAL(2.0, 2);
  TEST_EQUAL(2.0, 2.0f);
  TEST_EQUAL(2.0f, 2);
  TEST_EQUAL(2.0f, 2.0);
  enum Enum1 { A, B, C };
  enum Enum2 { AA, BB, CC };
  TEST_EQUAL(Enum1::A, Enum1::A);
  TEST_EQUAL(Enum1::A, Enum1::B); FAILURE_IS_SUCCESS;
  TEST_EQUAL(Enum1::A, Enum2::AA);
  TEST_EQUAL(Enum1::A, Enum2::BB); FAILURE_IS_SUCCESS;

  enum class Enum3 { A, B, C };
  enum class Enum4 { A, B, C };
  TEST_EQUAL(Enum3::A, Enum4::A);
  TEST_EQUAL(Enum3::A, Enum3::B); FAILURE_IS_SUCCESS;
  TEST_EQUAL(Enum3::B, Enum3::B);
  TEST_EQUAL(Enum3::C, Enum3::A); FAILURE_IS_SUCCESS;  
}
END_SECTION


START_SECTION("TEST_TRUE()")
{
  TEST_TRUE(2==3); FAILURE_IS_SUCCESS;
  TEST_TRUE(2==2);
}
END_SECTION

START_SECTION("TEST_FALSE()")
{
  TEST_FALSE(2 == 2);	FAILURE_IS_SUCCESS;
  TEST_FALSE(2 == 3);
}
END_SECTION

START_SECTION("TEST_REAL_SIMILAR()")
{

	const double b0 =  0.0;
	const double bn = -5.0;
	const double bp = +5.0;

	const double e0 =  0.0; // zero eps
	const double en = -0.1; // negative eps
	const double ep = +0.1; // positive eps
	double e;
	double f;

	std::string tmp_file_name;
	NEW_TMP_FILE(tmp_file_name);
	std::ofstream tmp_file(tmp_file_name.c_str());
	STATUS('\n' << tmp_file_name << ":0:  output of TEST_REAL_SIMILAR() elementary tests starts here");

#undef stdcout
#define stdcout tmp_file
	// The many {} are intended for code folding.  Do not mess them up.
	{
		{
			f = e0;
			{
				TOLERANCE_ABSOLUTE(0.0);
				TOLERANCE_RELATIVE(1.0);
				{
					{
						e = e0;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = ep;
						TEST_REAL_SIMILAR(b0+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f); FAILURE_IS_SUCCESS;
					}
					{
						e = en;
						TEST_REAL_SIMILAR(b0+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f); FAILURE_IS_SUCCESS;
					}
				}
			}
			{
				TOLERANCE_ABSOLUTE(0.25);
				TOLERANCE_RELATIVE(1.0);
				{
					{
						e = e0;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = ep;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = en;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
				}
			}
			{
				TOLERANCE_ABSOLUTE(0.0);
				TOLERANCE_RELATIVE(1.1);
				{
					{
						e = e0;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = ep;
						TEST_REAL_SIMILAR(b0+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = en;
						TEST_REAL_SIMILAR(b0+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
				}
			}
			{
				TOLERANCE_ABSOLUTE(0.25);
				TOLERANCE_RELATIVE(1.1);
// we need to get rid of some preprocessor macros, as the Microsoft VC++ compiler will crash with a stack overflow during
// compilation of this test otherwise - the stack size is not fixable apparently
#ifndef OPENMS_WINDOWSPLATFORM
				{
					{
						e = e0;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = ep;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = en;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
				}
#endif
			}
		}
		{
			f = en;
			{
				TOLERANCE_ABSOLUTE(0.0);
				TOLERANCE_RELATIVE(1.0);
				{
					{
						e = e0;
						TEST_REAL_SIMILAR(b0+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = ep;
						TEST_REAL_SIMILAR(b0+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f); FAILURE_IS_SUCCESS;
					}
					{
						e = en;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
				}
			}
			{
				TOLERANCE_ABSOLUTE(0.25);
				TOLERANCE_RELATIVE(1.0);
				{
					{
						e = e0;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = ep;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = en;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
				}
			}
			{
				TOLERANCE_ABSOLUTE(0.0);
				TOLERANCE_RELATIVE(1.1);
				{
					{
						e = e0;
						TEST_REAL_SIMILAR(b0+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = ep;
						TEST_REAL_SIMILAR(b0+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = en;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
				}
			}
			{
				TOLERANCE_ABSOLUTE(0.25);
				TOLERANCE_RELATIVE(1.1);
				{
					{
						e = e0;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = ep;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = en;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
				}
			}
		}
		{
			f = ep;
			{
				TOLERANCE_ABSOLUTE(0.0);
				TOLERANCE_RELATIVE(1.0);
				{
					{
						e = e0;
						TEST_REAL_SIMILAR(b0+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = ep;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = en;
						TEST_REAL_SIMILAR(b0+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f); FAILURE_IS_SUCCESS;
					}
				}
			}
			{
				TOLERANCE_ABSOLUTE(0.25);
				TOLERANCE_RELATIVE(1.0);
				{
					{
						e = e0;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = ep;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = en;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
				}
			}
			{
				TOLERANCE_ABSOLUTE(0.0);
				TOLERANCE_RELATIVE(1.1);
				{
					{
						e = e0;
						TEST_REAL_SIMILAR(b0+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = ep;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = en;
						TEST_REAL_SIMILAR(b0+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
				}
			}
			{
				TOLERANCE_ABSOLUTE(0.25);
				TOLERANCE_RELATIVE(1.1);
				{
					{
						e = e0;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = ep;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
					{
						e = en;
						TEST_REAL_SIMILAR(b0+e,b0+f);
						TEST_REAL_SIMILAR(b0+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(b0+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bn+e,bn+f);
						TEST_REAL_SIMILAR(bn+e,bp+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,b0+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bn+f); FAILURE_IS_SUCCESS;
						TEST_REAL_SIMILAR(bp+e,bp+f);
					}
				}
			}
		}
	}
#undef stdcout
#define stdcout std::cout

}
END_SECTION

#if 1

START_SECTION("TEST_STRING_SIMILAR")
{

	const char lhs[] = "a  bcd  ef 10.0 ghi jk   l\n l 101.125mno p \nqrs";
	const char rhs[] = "a \t bcd ef 12.0 ghi  jk l\n l 124.125mno  p  \nqrs";
	TOLERANCE_ABSOLUTE(1.);
	TOLERANCE_RELATIVE(1.3);
	TEST_STRING_SIMILAR(lhs,rhs);
	TOLERANCE_ABSOLUTE(30.);
	TOLERANCE_RELATIVE(1.1);
	TEST_STRING_SIMILAR(lhs,rhs);

	//--------------------------------------------------------

	const double numbers[] = { -5.1, -5.0, -4.9, -0.1, 0.0, 0.1, 4.9, 5.0, 5.1 };
	UInt number_of_numbers = sizeof(numbers)/sizeof(*numbers);
	std::vector<OpenMS::String> number_strings;
	for ( UInt i = 0; i < number_of_numbers; ++i )
	{
		number_strings.push_back(OpenMS::String("ABC") + numbers[i] + "XYZ");
	}

	const double tolerance_absolute[2] = { 0.0, 0.25 };
	const double tolerance_relative[2] = { 1.0, 1.1 };

	// Debugging.  If you really want to know it.  Output is > 10000 lines.
	const bool compare_always = false;

	for ( UInt ta = 0; ta < 2; ++ta )
	{
		TOLERANCE_ABSOLUTE(tolerance_absolute[ta]);
		for ( UInt tr = 0; tr < 2; ++tr )
		{
			TOLERANCE_RELATIVE(tolerance_relative[tr]);

			for ( UInt i = 0; i < number_of_numbers; ++i )
			{
				const double ni = numbers[i];
				const OpenMS::String& si = number_strings[i];

				for ( UInt j = 0; j < number_of_numbers; ++j )
				{
					const double nj = numbers[j];
					const OpenMS::String& sj = number_strings[j];

          // Bypass the macros to avoid lengthy output. These functions do the real job.
          bool save = TEST::test;
          const bool ne = TEST::isRealSimilar(ni,nj);
          TEST::testStringSimilar(__FILE__,__LINE__,si,"si",sj,"sj");
          const bool se = TEST::this_test;
          TEST::this_test = true;
          TEST::test = save;

					if ( se != ne || compare_always )
					{
						// We have an issue. Get the message.
						STATUS(" ni:" << ni << "  nj:" << nj << "  si:" << si << "  sj:" << sj);
						 // The real question.
						TEST_EQUAL(se,ne);

						// The next two TEST_.. should produce the same decision and similar messages.
						bool save = TEST::test;
						TEST_REAL_SIMILAR(ni,nj); // should be equal to ne
						TEST_STRING_SIMILAR(si,sj); // should be equal to se
						TEST::test = save;
					}

				}
			}
		}
	}

}
END_SECTION
#endif

#if 1

START_SECTION("TEST_FILE_SIMILAR")
{

	std::string filename1, filename2;
	NEW_TMP_FILE(filename1);
	NEW_TMP_FILE(filename2);
	{
		std::ofstream file1(filename1.c_str());
		std::ofstream file2(filename2.c_str());
		file1 << "1 \n xx\n 2.008	\n 3" << std::flush;
		file2 << "1.08 \n    xx\n		\n\n  				  	0002.04000 \n 3" << std::flush;
		file1.close();
		file2.close();
	}

	TOLERANCE_ABSOLUTE(0.01);
	TOLERANCE_RELATIVE(1.1);
	TEST_FILE_SIMILAR(filename1,filename2);
}
END_SECTION

START_SECTION("TEST_EQUAL")
	TEST_EQUAL(1.0, 1.0)
	TEST_EQUAL('A', 'A')
END_SECTION

START_SECTION("TEST_NOT_EQUAL")
	TEST_NOT_EQUAL(0, 1)
	TEST_NOT_EQUAL('A', 'B')
END_SECTION

START_SECTION("TEST_EXCEPTION")
	TEST_EXCEPTION(Exception::NullPointer, throw Exception::NullPointer(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION))
END_SECTION

START_SECTION("TEST_EXCEPTION_WITH_MESSAGE")
	TEST_EXCEPTION_WITH_MESSAGE(Exception::NullPointer, throw Exception::NullPointer(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION), "a null pointer was specified")
END_SECTION

START_SECTION("TEST_PRECONDITION_VIOLATED")
  // recommended usage, success
  TEST_PRECONDITION_VIOLATED(throw_a_Precondition_Exception());
  int this_was_evaluated = false;
  // recommended usage, but failure will be signaled only when compiled in Debug mode.
  TEST_PRECONDITION_VIOLATED(this_was_evaluated = true);  if ( this_was_evaluated )  { FAILURE_IS_SUCCESS; }
  // wrong exception thrown, or none at all
  TEST_PRECONDITION_VIOLATED(throw_a_Postcondition_Exception()); if ( this_was_evaluated )  { FAILURE_IS_SUCCESS; }

#ifndef OPENMS_ASSERTIONS
  NOT_TESTABLE; // just to avoid a warning message in Release mode - all test macros will expand empty.
#endif

  END_SECTION

START_SECTION("TEST_POSTCONDITION_VIOLATED")
  // recommended usage, success
  TEST_POSTCONDITION_VIOLATED(throw_a_Postcondition_Exception());
  int this_was_evaluated = false;
  // recommended usage, but failure will be signaled only when compiled in Debug mode.
  TEST_POSTCONDITION_VIOLATED(this_was_evaluated = true);  if ( this_was_evaluated )  { FAILURE_IS_SUCCESS; }
  // wrong exception thrown, or none at all
  TEST_POSTCONDITION_VIOLATED(throw_a_Precondition_Exception()); if ( this_was_evaluated )  { FAILURE_IS_SUCCESS; }

#ifndef OPENMS_ASSERTIONS
  NOT_TESTABLE; // just to avoid a warning message in Release mode - all test macros will expand empty.
#endif

  END_SECTION

START_SECTION("OPENMS_PRETTY_FUNCTION")
	struct Dummy
	{
		std::string f_dummy(double, float,int,unsigned,long,unsigned long,char) { return OPENMS_PRETTY_FUNCTION; }
	} dummy;
	STATUS("\n\n\tExample for usage of OPENMS_PRETTY_FUNCTION inside a member function of a nested class in main():\n\t" << dummy.f_dummy(0,0,0,0,0,0,0) << '\n')
	NOT_TESTABLE
END_SECTION

START_SECTION("STATUS")
	STATUS("status message")
	NOT_TESTABLE
END_SECTION

START_SECTION("TEST_FILE_EQUAL")
	TEST_FILE_EQUAL(OPENMS_GET_TEST_DATA_PATH("class_test_infile.txt"), OPENMS_GET_TEST_DATA_PATH("class_test_template.txt"))
END_SECTION

START_SECTION("ABORT_IF")
	TEST_EQUAL(1, 1)
	while (true)
  {
    ABORT_IF(true) // will 'break;' internally, but we do not want to leave the test
	}
  FAILURE_IS_SUCCESS;			 	
	
END_SECTION

START_SECTION("TEST_REAL_SIMILAR : type checking")
{
 	TEST_REAL_SIMILAR( 0.0  , 0.0  );
 	TEST_REAL_SIMILAR( 0.0  , 0.0F );
 	TEST_REAL_SIMILAR( 0.0  , 0.0L );
 	TEST_REAL_SIMILAR( 0.0F , 0.0  );
 	TEST_REAL_SIMILAR( 0.0F , 0.0F );
 	TEST_REAL_SIMILAR( 0.0F , 0.0L );
 	TEST_REAL_SIMILAR( 0.0L , 0.0  );
 	TEST_REAL_SIMILAR( 0.0L , 0.0F );
 	TEST_REAL_SIMILAR( 0.0L , 0.0L );

	TEST_REAL_SIMILAR( 0.0  , 0U   );
	TEST_REAL_SIMILAR( 0U   , 0.0  ); FAILURE_IS_SUCCESS;
	TEST_REAL_SIMILAR( 0U   , 0U   ); FAILURE_IS_SUCCESS;

	TEST_REAL_SIMILAR( 0.0  , 0L   );
	TEST_REAL_SIMILAR( 0L   , 0.0  ); FAILURE_IS_SUCCESS;
	TEST_REAL_SIMILAR( 0L   , 0L   ); FAILURE_IS_SUCCESS;

	TEST_REAL_SIMILAR( 0   , 0U   ); FAILURE_IS_SUCCESS;
	TEST_REAL_SIMILAR( 0   , 0L   ); FAILURE_IS_SUCCESS;
	TEST_REAL_SIMILAR( 0   , 0UL  ); FAILURE_IS_SUCCESS;

	TEST_REAL_SIMILAR( 0.0  , 0UL  );
	TEST_REAL_SIMILAR( 0UL  , 0.0  ); FAILURE_IS_SUCCESS;
	TEST_REAL_SIMILAR( 0UL  , 0UL  ); FAILURE_IS_SUCCESS;

	TEST_REAL_SIMILAR( 0.0F , 0    );
	TEST_REAL_SIMILAR( 0.0  , 0    );
	TEST_REAL_SIMILAR( 0.0L , 0    );
	TEST_REAL_SIMILAR( 0    , 0.0F ); FAILURE_IS_SUCCESS;
	TEST_REAL_SIMILAR( 0    , 0.0  ); FAILURE_IS_SUCCESS;
	TEST_REAL_SIMILAR( 0    , 0.0L ); FAILURE_IS_SUCCESS;

	TEST_REAL_SIMILAR( 0.0F , 0U   );
	TEST_REAL_SIMILAR( 0.0  , 0U   );
	TEST_REAL_SIMILAR( 0.0L , 0U   );
	TEST_REAL_SIMILAR( 0U   , 0.0F ); FAILURE_IS_SUCCESS;
	TEST_REAL_SIMILAR( 0U   , 0.0  ); FAILURE_IS_SUCCESS;
	TEST_REAL_SIMILAR( 0U   , 0.0L ); FAILURE_IS_SUCCESS;

	TEST_REAL_SIMILAR( std::numeric_limits<double>::quiet_NaN(), 0.0 ); FAILURE_IS_SUCCESS;
  TEST_REAL_SIMILAR( 0.0, std::numeric_limits<double>::quiet_NaN() ); FAILURE_IS_SUCCESS;
  TEST_REAL_SIMILAR( std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN() ); FAILURE_IS_SUCCESS;

  TEST_REAL_SIMILAR( std::numeric_limits<float>::quiet_NaN(), 0.0 ); FAILURE_IS_SUCCESS;
  TEST_REAL_SIMILAR( 0.0, std::numeric_limits<float>::quiet_NaN() ); FAILURE_IS_SUCCESS;
  TEST_REAL_SIMILAR( std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN() ); FAILURE_IS_SUCCESS;

  TEST_REAL_SIMILAR( std::numeric_limits<long double>::quiet_NaN(), 0.0 ); FAILURE_IS_SUCCESS;
  TEST_REAL_SIMILAR( 0.0, std::numeric_limits<long double>::quiet_NaN() ); FAILURE_IS_SUCCESS;
  TEST_REAL_SIMILAR( std::numeric_limits<long double>::quiet_NaN(), std::numeric_limits<long double>::quiet_NaN() ); FAILURE_IS_SUCCESS;

}
END_SECTION

#endif

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
