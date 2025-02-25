// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

#include <OpenMS/CONCEPT/FuzzyStringComparator.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <fstream>

/////////////////////////////////////////////////////////////

START_TEST(FuzzyStringComparator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FuzzyStringComparator* inst_ptr = nullptr;
FuzzyStringComparator* nullPointer = nullptr;
START_SECTION((FuzzyStringComparator()))
{
	inst_ptr = new FuzzyStringComparator;
  TEST_NOT_EQUAL(inst_ptr, nullPointer);
}
END_SECTION

START_SECTION((virtual ~FuzzyStringComparator()))
{
	delete inst_ptr;
}
END_SECTION

START_SECTION((FuzzyStringComparator& operator=(const FuzzyStringComparator& rhs)))
{
  // Not implemented
	NOT_TESTABLE;
}
END_SECTION

START_SECTION((FuzzyStringComparator(const FuzzyStringComparator& rhs)))
{
  // Not implemented
	NOT_TESTABLE;
}
END_SECTION

//------------------------------------------------------------

START_SECTION((const double& getAcceptableAbsolute() const))
{
	// tested along with set-method
	NOT_TESTABLE;
}
END_SECTION

START_SECTION((const double& getAcceptableRelative() const))
{
	// tested along with set-method
	NOT_TESTABLE;
}
END_SECTION

START_SECTION((const int& getVerboseLevel() const))
{
	// tested along with set-method
	NOT_TESTABLE;
}
END_SECTION

START_SECTION((const int& getTabWidth() const))
{
	// tested along with set-method
	NOT_TESTABLE;
}
END_SECTION

START_SECTION((const int& getFirstColumn() const))
{
	// tested along with set-method
	NOT_TESTABLE;
}
END_SECTION

START_SECTION((std::ostream& getLogDestination() const))
{
	// tested along with set-method
	NOT_TESTABLE;
}
END_SECTION

START_SECTION((void setAcceptableAbsolute(const double rhs)))
{
	FuzzyStringComparator fsc;
	fsc.setAcceptableAbsolute(2345.6789);
	TEST_REAL_SIMILAR(fsc.getAcceptableAbsolute(),2345.6789);
}
END_SECTION

START_SECTION((void setAcceptableRelative(const double rhs)))
{
	FuzzyStringComparator fsc;
	fsc.setAcceptableRelative(6789.2345);
	TEST_REAL_SIMILAR(fsc.getAcceptableRelative(),6789.2345);
}
END_SECTION

START_SECTION((void setTabWidth(const int rhs)))
{
	FuzzyStringComparator fsc;
	fsc.setTabWidth(1452);
	TEST_EQUAL(fsc.getTabWidth(),1452);
}
END_SECTION

START_SECTION((void setFirstColumn(const int rhs)))
{
	FuzzyStringComparator fsc;
	fsc.setFirstColumn(4321235);
	TEST_EQUAL(fsc.getFirstColumn(),4321235);
}
END_SECTION

START_SECTION((void setLogDestination(std::ostream & rhs)))
{
	FuzzyStringComparator fsc;
	TEST_EQUAL(&fsc.getLogDestination(),&std::cout);
	fsc.setLogDestination(std::cerr);
	TEST_EQUAL(&fsc.getLogDestination(),&std::cerr);
	TEST_NOT_EQUAL(&fsc.getLogDestination(),&std::cout);
	fsc.setLogDestination(std::cout);
	TEST_NOT_EQUAL(&fsc.getLogDestination(),&std::cerr);
	TEST_EQUAL(&fsc.getLogDestination(),&std::cout);
}
END_SECTION

START_SECTION((void setVerboseLevel(const int rhs)))
{
	FuzzyStringComparator fsc;
	// default should be 2
	TEST_EQUAL(fsc.getVerboseLevel(),2);
	fsc.setVerboseLevel(88);
	TEST_EQUAL(fsc.getVerboseLevel(),88);
	fsc.setVerboseLevel(-21);
	TEST_EQUAL(fsc.getVerboseLevel(),-21);
}
END_SECTION

{
  FuzzyStringComparator fsc;
  FuzzyStringComparator const & fsc_cref = fsc;

  START_SECTION((const StringList& getWhitelist() const ))
  {
    TEST_EQUAL(fsc_cref.getWhitelist().empty(),true);
    // continued below
  }
  END_SECTION

  START_SECTION((StringList& getWhitelist()))
  {
    TEST_EQUAL(fsc_cref.getWhitelist().empty(),true);
    // continued below
  }
  END_SECTION

  START_SECTION((void setWhitelist(const StringList &rhs)))
  {
    fsc.setWhitelist(ListUtils::create<String>("null,eins,zwei,drei"));
    TEST_STRING_EQUAL(fsc.getWhitelist()[0],"null");
    TEST_STRING_EQUAL(fsc_cref.getWhitelist()[1],"eins");
    TEST_EQUAL(fsc_cref.getWhitelist().size(),4);
    fsc.setWhitelist(ListUtils::create<String>("zero,one,two,three,four"));
    TEST_STRING_EQUAL(fsc.getWhitelist()[0],"zero");
    TEST_STRING_EQUAL(fsc_cref.getWhitelist()[1],"one");
    TEST_EQUAL(fsc_cref.getWhitelist().size(),5);
  }
  END_SECTION

}


//------------------------------------------------------------

START_SECTION((bool compareStrings(std::string const &lhs, std::string const &rhs)))
{
	std::ostringstream log;
	//------------------------------
	// A few test to show what regular expressions could not do but our class can do.
	{
		FuzzyStringComparator fsc;
		log.str("");
		fsc.setLogDestination(log);
		fsc.setVerboseLevel(2);
		fsc.setAcceptableRelative(1.00021);
		fsc.setAcceptableAbsolute(0.0);
		bool result = (fsc.compareStrings("0.9999E4","1.0001E4")!=0);
		TEST_EQUAL(result,true);
		// STATUS(log.str());
	}
	{
		FuzzyStringComparator fsc;
		log.str("");
		fsc.setLogDestination(log);
		fsc.setVerboseLevel(2);
		fsc.setAcceptableRelative(1.0);
		fsc.setAcceptableAbsolute(2.0);
		bool result = (fsc.compareStrings("0.9999E4","1.0001E4")!=0);
		TEST_EQUAL(result,true);
		// STATUS(log.str());
	}
	//------------------------------
	// Various issues, mixing letters, whitespace, and numbers.
	{
		FuzzyStringComparator fsc;
		log.str("");
		fsc.setLogDestination(log);
		fsc.setVerboseLevel(1);
		fsc.setAcceptableRelative(1.01);
		fsc.setAcceptableAbsolute(0.001);
		bool result = (fsc.compareStrings("bl   a b 00.0022 asdfdf","bl a  b 0.00225 asdfdf")!=0);
		TEST_EQUAL(result,true);
		// STATUS(log.str());
	}
	{
		FuzzyStringComparator fsc;
		log.str("");
		fsc.setLogDestination(log);
		fsc.setVerboseLevel(1);
		fsc.setAcceptableRelative(1.01);
		fsc.setAcceptableAbsolute(0.01);
		bool result = (fsc.compareStrings("bl   a 1.2   b","bl a 1.25 b")!=0);
		TEST_EQUAL(result,false);
		// STATUS(log.str());
	}
	{
		FuzzyStringComparator fsc;
		log.str("");
		fsc.setLogDestination(log);
		fsc.setVerboseLevel(1);
		fsc.setAcceptableRelative(2.);
		fsc.setAcceptableAbsolute(0.01);
		bool result = (fsc.compareStrings("bl   a 1.2   b","bl a 1.25 b")!=0);
		TEST_EQUAL(result,true);
		// STATUS(log.str());
	}
	{
		FuzzyStringComparator fsc;
		log.str("");
		fsc.setLogDestination(log);
		fsc.setVerboseLevel(1);
		fsc.setAcceptableRelative(1.01);
		fsc.setAcceptableAbsolute(0.0);
		bool result = (fsc.compareStrings("bl   a 1.002   b","bl a 1.0025 b")!=0);
		TEST_EQUAL(result,true);
		// STATUS(log.str());
	}
	//------------------------------
	// Test the impact of verbosity_level.
	{
		FuzzyStringComparator fsc;
		log.str("");
		fsc.setLogDestination(log);
		fsc.setVerboseLevel(1);
		fsc.setAcceptableRelative(1.03);
		fsc.setAcceptableAbsolute(0.01);
		bool result = (fsc.compareStrings("1 \n 		   2	\n 3","1.01 \n \n		\n\n  					  	0002.01000 \n 3")!=0);
		TEST_EQUAL(result,true);
		std::vector<OpenMS::String> substrings;
		result = OpenMS::String(log.str()).split('\n',substrings);
		// STATUS(log.str());
		TEST_EQUAL(result, false);
		TEST_EQUAL(substrings.size(), 0);
		STATUS(substrings.size());
	}
	{
		FuzzyStringComparator fsc;
		log.str("");
		fsc.setLogDestination(log);
		fsc.setVerboseLevel(2);
		fsc.setAcceptableRelative(1.03);
		fsc.setAcceptableAbsolute(0.01);
		bool result = (fsc.compareStrings("1 \n 		   2	\n 3","1.01 \n \n		\n\n  					  	0002.01000 \n 3")!=0);
		TEST_EQUAL(result,true);
		std::vector<OpenMS::String> substrings;
		OpenMS::String(log.str()).split('\n',substrings);
		// STATUS(log.str());
		// Magic alert! - You might need to edit these numbers if reportSuccess_() or reportFailure_() changes.
		TEST_EQUAL(substrings.size(),17);
		ABORT_IF(substrings.size()!=17);
		TEST_STRING_EQUAL(substrings[0],"PASSED.");
	}
	{
		FuzzyStringComparator fsc;
		log.str("");
		fsc.setLogDestination(log);
		fsc.setVerboseLevel(1);
		fsc.setAcceptableRelative(1.01);
		fsc.setAcceptableAbsolute(0.01);
		fsc.compareStrings("1 \n 		   2	\n 3","1.11 \n \n		\n\n  					  	0004.01000 \n 3");
		std::vector<OpenMS::String> substrings;
		OpenMS::String(log.str()).split('\n',substrings);
		// STATUS(log.str());
		// Magic alert! - You might need to edit these numbers if reportSuccess_() or reportFailure_() changes.
		TEST_EQUAL(substrings.size(),36);
		ABORT_IF(substrings.size()!=36);
		TEST_STRING_EQUAL(substrings[0],"FAILED: 'ratio of numbers is too large'");
	}
	{
		FuzzyStringComparator fsc;
		log.str("");
		fsc.setLogDestination(log);
		//		fsc.setLogDestination(std::cout);
		fsc.setVerboseLevel(3);
		fsc.setAcceptableRelative(1.01);
		fsc.setAcceptableAbsolute(0.01);
		fsc.compareStrings
			(
			 "1 \n xx\n 2.008	\n 3",
			 "1.11 \nU\n		\n\n  q					  	0002.04000 \n 3"
			);
		std::vector<OpenMS::String> substrings;
		OpenMS::String(log.str()).split('\n',substrings);
		// STATUS(log.str());
		// Magic alert! - You might need to edit these numbers if reportSuccess_() or reportFailure_() changes.
		TEST_EQUAL(substrings.size(),246);
		ABORT_IF(substrings.size()!=246);
		TEST_STRING_EQUAL(substrings[0],"FAILED: 'ratio of numbers is too large'");
		TEST_STRING_EQUAL(substrings[35],"FAILED: 'input_1 is whitespace, but input_2 is not'");
		TEST_STRING_EQUAL(substrings[70],"FAILED: 'different letters'");
		TEST_STRING_EQUAL(substrings[105],"FAILED: 'line from input_2 is shorter than line from input_1'");
		TEST_STRING_EQUAL(substrings[140],"FAILED: 'input_1 is a number, but input_2 is not'");
		TEST_STRING_EQUAL(substrings[175],"FAILED: 'input_1 is not a number, but input_2 is'");
		TEST_STRING_EQUAL(substrings[210],"FAILED: 'line from input_1 is shorter than line from input_2'");
	}
  {
		FuzzyStringComparator fsc;
		log.str("");
		fsc.setLogDestination(log);
		fsc.setVerboseLevel(2);
		fsc.setAcceptableRelative(1.0);
		fsc.setAcceptableAbsolute(2.0);
		bool result = fsc.compareStrings("0.9999X","1.0001X");
		TEST_EQUAL(result,true);
		// STATUS(log.str());
  }
}
END_SECTION


START_SECTION((bool compareStreams(std::istream &input_1, std::istream &input_2)))
{
	std::ostringstream log;
	{
		FuzzyStringComparator fsc;
		log.str("");
		fsc.setLogDestination(log);
		fsc.setVerboseLevel(3);
		fsc.setAcceptableRelative(1.01);
		fsc.setAcceptableAbsolute(0.01);
		std::istringstream lhs("1 \n xx\n 2.008	\n 3");
		std::istringstream rhs("1.11 \nU\n		\n\n  q					  	0002.04000 \n 3");
		fsc.compareStreams(lhs,rhs);
		std::vector<OpenMS::String> substrings;
		OpenMS::String(log.str()).split('\n',substrings);
		// STATUS(log.str());
		// Magic alert! - You might need to edit these numbers if reportSuccess_() or reportFailure_() changes.
		TEST_EQUAL(substrings.size(),246);
		ABORT_IF(substrings.size()!=246);
		TEST_STRING_EQUAL(substrings[0],"FAILED: 'ratio of numbers is too large'");
		TEST_STRING_EQUAL(substrings[35],"FAILED: 'input_1 is whitespace, but input_2 is not'");
		TEST_STRING_EQUAL(substrings[70],"FAILED: 'different letters'");
		TEST_STRING_EQUAL(substrings[105],"FAILED: 'line from input_2 is shorter than line from input_1'");
		TEST_STRING_EQUAL(substrings[140],"FAILED: 'input_1 is a number, but input_2 is not'");
		TEST_STRING_EQUAL(substrings[175],"FAILED: 'input_1 is not a number, but input_2 is'");
		TEST_STRING_EQUAL(substrings[210],"FAILED: 'line from input_1 is shorter than line from input_2'");
	}
}
END_SECTION

START_SECTION((bool compareFiles(const std::string &filename_1, const std::string &filename_2)))
{
	std::ostringstream log;
	{
		FuzzyStringComparator fsc;
		log.str("");
		fsc.setLogDestination(log);
		fsc.setVerboseLevel(3);
		fsc.setAcceptableRelative(1.01);
		fsc.setAcceptableAbsolute(0.01);
		std::string filename1, filename2;
		NEW_TMP_FILE(filename1);
		NEW_TMP_FILE(filename2);
		{
			std::ofstream file1(filename1.c_str());
			std::ofstream file2(filename2.c_str());
			file1 << "1 \n xx\n 2.008	\n 3" << std::flush;
			file2 << "1.11 \nU\n		\n\n  q					  	0002.04000 \n 3" << std::flush;
			file1.close();
			file2.close();
		}
		fsc.compareFiles(filename1,filename2);
		std::vector<OpenMS::String> substrings;
		OpenMS::String(log.str()).split('\n',substrings);
		// STATUS(log.str());
		// Magic alert! - You might need to edit these numbers if reportSuccess_() or reportFailure_() changes.
		TEST_EQUAL(substrings.size(),246);
		ABORT_IF(substrings.size()!=246);
		TEST_STRING_EQUAL(substrings[0],"FAILED: 'ratio of numbers is too large'");
		TEST_STRING_EQUAL(substrings[35],"FAILED: 'input_1 is whitespace, but input_2 is not'");
		TEST_STRING_EQUAL(substrings[70],"FAILED: 'different letters'");
		TEST_STRING_EQUAL(substrings[105],"FAILED: 'line from input_2 is shorter than line from input_1'");
		TEST_STRING_EQUAL(substrings[140],"FAILED: 'input_1 is a number, but input_2 is not'");
		TEST_STRING_EQUAL(substrings[175],"FAILED: 'input_1 is not a number, but input_2 is'");
		TEST_STRING_EQUAL(substrings[210],"FAILED: 'line from input_1 is shorter than line from input_2'");
	}
}
END_SECTION

//------------------------------------------------------------

// START_SECTION(void reportFailure_( char const * const message ) const throw(Failure))
// {
// 	// Tested in compare...() methods
//   NOT_TESTABLE;
// }
// END_SECTION

// START_SECTION(void reportSuccess_() const)
// {
// 	// Tested in compare...() methods
//   NOT_TESTABLE;
// }
// END_SECTION

END_TEST
