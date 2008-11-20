// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/DATASTRUCTURES/String.h>

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
			std__cout << __FILE__ ":" <<  __LINE__ <<													\
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
			std__cout << __FILE__ ":" <<  __LINE__ <<													\
				" error:  The preceeding test was supposed to fail, but it did not.  =>  FAILURE" << \
				std::endl;																											\
		}																																		\
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

START_SECTION()
	std::string tmp_filename;
	NEW_TMP_FILE(tmp_filename);
	TEST::this_test = (tmp_filename != "");
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

#undef std__cout
#define std__cout tmp_file
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
#undef std__cout
#define std__cout std::cout

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
					const bool ne = TEST::isRealSimilar(ni,nj);
					const bool se = TEST::isStringSimilar(si,sj);

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

START_SECTION()
	TEST_EQUAL(1.0, 1.0)
	TEST_EQUAL('A', 'A')
END_SECTION

START_SECTION()
	TEST_NOT_EQUAL(0, 1)
	TEST_NOT_EQUAL('A', 'B')
END_SECTION

START_SECTION()
	TEST_EXCEPTION(Exception::NullPointer, throw Exception::NullPointer(__FILE__, __LINE__, __PRETTY_FUNCTION__))
END_SECTION

START_SECTION()
	TEST_EXCEPTION_WITH_MESSAGE(Exception::NullPointer, throw Exception::NullPointer(__FILE__, __LINE__, __PRETTY_FUNCTION__), "a null pointer was specified")
END_SECTION

START_SECTION()
	struct Dummy
	{
		std::string f_dummy(double, float,int,unsigned,long,unsigned long,char) { return __PRETTY_FUNCTION__; }
	} dummy;
	STATUS("\n\n\tExample for usage of __PRETTY_FUNCTION__ inside a member function of a nested class in main():\n\t" << dummy.f_dummy(0,0,0,0,0,0,0) << '\n')
END_SECTION

START_SECTION()
	STATUS("status message")
END_SECTION

START_SECTION()
	TEST_FILE_EQUAL("data/class_test_infile.txt", "data/class_test_template.txt")
END_SECTION

START_SECTION()
	ABORT_IF(true)
	TEST_EQUAL(1, 0)
END_SECTION

#endif

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
