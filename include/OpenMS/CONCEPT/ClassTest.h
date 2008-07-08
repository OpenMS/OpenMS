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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_CLASSTEST_H
#define OPENMS_CONCEPT_CLASSTEST_H

/** @brief Indicates that a class test is being compiled.

	Used e.g. in OPENMS_PRECONDITION and OPENMS_POSTCONDITION so that	we
	can	test these even if the global OPENMS_DEBUG macro is not set.
*/
#define OPENMS_WITHIN_CLASSTEST 1

// Avoid OpenMS includes here at all costs
// When the included headers are changed, *all* tests have to be recompiled!
// Use the ClassTest class if you need add high-level functionality.
// Includes in the C-file are ok...
#include <OpenMS/config.h>
#include <OpenMS/SYSTEM/ProcessResource.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <vector>
#include <string>
#include <cstring>
#include <list>
#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>

#include <unistd.h> // unlink()
#include <stdio.h>  // tmpnam()
#include <math.h>   // fabs
#include <stdlib.h> // getenv()

#ifdef OPENMS_HAS_SSTREAM
# include <sstream>
#else
# include <strstream>
#endif

// Empty declaration to avoid problems in case the namespace is not
// yet defined (e.g. TEST/ClassTest_test.C)
namespace OpenMS
{
	namespace Internal
	{
		/// Auxilary class for class tests
		class ClassTest
		{
			
			public:
			
				/**
					@brief Valides the given files against the XML schema (if available)
					@return If all files passed the validation 
				*/
				static bool validate(const std::vector<std::string>& file_names);
				
				///Creates a temporary file name from the test name and the line
				static std::string tmpFileName(const std::string& file, int line);
		};
	}
}

/**
	@defgroup ClassTest Class test macros

  @brief Macros used by the test programs in the subdirectory <code>OpenMS/source/TEST</code>.

	On successful operation the test program will print out the message "PASSED",
	otherwise "FAILED".

	If called with the -v option, the test program prints verbose information
	about individual tests.

	If called with the -V option, the test program prints even more verbose
	information for every subtest.

	@ingroup Concept

	@internal Emacs has a fantastic command "backslashify" that you can use to
	line up the backslashes in the define blocks.
*/
//@{

/**	@brief Define the precision for floating point comparisons.

  The macro #TEST_REAL_EQUAL checks whether the floating point number returned
  by the subtest is close to the expected result by comparing the absolute
  value of the difference of the two values to #PRECISION.

	The default value is $10^{-6}$. It is possible to redefine precision in the
	test program by calling this macro with the new value.

	 @hideinitializer
*/
#define PRECISION(a) \
		TEST::precision = (a);																														\
    if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))							\
    {																																									\
			if (!TEST::newline)																															\
			{																																								\
				TEST::newline = true;																													\
				std::cout << std::endl;																												\
			}																																								\
      std::cout << "    (line " << __LINE__ << ":  PRECISION(" << TEST::precision			\
								<< ")   ("#a")" << std::endl; 																				\
    }																																									\


/**	@brief Create the test header for a certain class.

  This macro defines the start of the test program for a given classname.  The
  classname is printed together with some information when calling the test
  program with any arguments (except for <code>-v</code> or <code>-V</code>).

	This macro should be the first to call in a test program. It introduces a
	global <code>try</code> block to catch any unwanted exceptions. If any of
	these exceptions occurs, all tests failed.  Exceptions defined by OpenMS
	(i.e. exception classes derived from Exception::BaseException) provide some
	additional information that is evaluated by the #END_TEST macro. The
	#END_TEST macro also closes the <code>try</code> block. This <code>try</code>
	block should never catch an exception! All exceptions that are thrown due to
	some malfunction in one of the member functions should be caught by the
	<code>try</code> block created by #CHECK and #RESULT.

   @hideinitializer
*/
#define START_TEST(class_name, version)																							\
/* define a special namespace for all internal variables */													\
/* to avoid potential collisions                         */													\
namespace TEST {																																		\
	int											verbose = 0;																							\
	bool										all_tests = true;																					\
  bool										test = true;																							\
	bool										this_test;																								\
	int											exception = 0;																						\
	std::string							exception_name = "";																			\
	std::string							exception_message = "";																		\
	std::string							test_name = "";																						\
  int											check_line = 0;																						\
  int											test_line = 0;																						\
	const char*							version_string = version;																	\
	bool										newline = false;																					\
	std::vector<std::string>	tmp_file_list;																					\
	std::ifstream						infile;																										\
	std::ifstream						templatefile;																							\
	bool										equal_files;																							\
	double									precision = 1e-5;																					\
	char										line_buffer[65537];                                       \
	int                     test_count = 0;                                           \
	std::string             add_message = "";  																		    \
}																																										\
																																										\
																																										\
int main(int argc, char **argv)																											\
{																																										\
																																										\
	if (argc == 2) {																																	\
		if (!strcmp(argv[1], "-v"))																											\
			TEST::verbose = 1;																														\
		if (!strcmp(argv[1], "-V"))																											\
			TEST::verbose = 2;																														\
	};																																								\
																																										\
	if ((argc > 2) || ((argc == 2) && (TEST::verbose == 0))) {												\
		std::cerr																																				\
     << "This is " << argv[0] << ", the test program for the " << std::endl					\
     << #class_name " class." << std::endl << std::endl															\
     << "On successful operation it simply returns PASSED," << std::endl						\
		 << "otherwise FAILED is printed." << std::endl																	\
		 << "If called with an argument of -v, " << std::endl														\
		 << "prints detailed information about individual tests." << std::endl					\
		 << "Option -V provides verbose information on every subtest." << std::endl;		\
		return 1;																																				\
	}																																									\
																																										\
	if (TEST::verbose > 0)																														\
		std::cout << "Version: " << TEST::version_string << std::endl;									\
																																										\
  char * pPath;                                                                     \
  pPath = getenv ("OPENMS_TESTTIMEOUT");                                            \
  int timeout = 300;                                                                 \
                                                                                    \
  if (pPath!=NULL)                                                                  \
  {                                                                                 \
    try                                                                             \
    {                                                                               \
      timeout = boost::lexical_cast<int>(pPath);                                    \
    }                                                                               \
    catch (boost::bad_lexical_cast&)                                                \
    {                                                                               \
      /* std::cout << "cast failed. timeout is: " << timeout << "\n\n";  */         \
    }                                                                               \
  }                                                                                 \
  OpenMS::ProcessResource::LimitCPUTime(timeout);                                   \
                                                                                    \
  try {


/**	@brief Termination of test program.

  This macro implements the correct termination of the test program and should
  therefore be the last macro to call.  It determines the exit code based on
  all previously run subtests and prints out the message "PASSED" or "FAILED".
  This macro also closes the global <code>try</code> block opened by
  #START_TEST and contains the related <code>catch</code> clauses. If an
  exception is caught here, the test program fails.

   @hideinitializer
*/
#define END_TEST																																		\
	/* global try block */																														\
	}																																									\
	/* catch FileNotFound exceptions to print out the file name */										\
	catch (OpenMS::Exception::FileNotFound e)																					\
	{																																									\
		TEST::this_test = false;																												\
		TEST::test = false;																															\
		TEST::all_tests = false;																												\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))						\
		{																																								\
			if (TEST::exception == 1) /* dummy to avoid compiler warnings */							\
				TEST::exception++;																													\
			std::cout << std::endl << "    (caught exception of type `"										\
			          << e.getName() << "'";																							\
			if ((e.getLine() > 0) && (std::strcmp(e.getFile(),"")!=0))															\
				std::cout << " outside a subtest, which was thrown in line " << e.getLine()	\
									<< " of file " << e.getFile()																			\
									<< " in function " << e.getFunction();														\
			std::cout << " - unexpected!) " << std::endl;																	\
		}																																								\
  }																																									\
	/* catch OpenMS exceptions to retrieve additional information */									\
	catch (OpenMS::Exception::BaseException& e)																								\
	{																																									\
		TEST::this_test = false;																												\
		TEST::test = false;																															\
		TEST::all_tests = false;																												\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))						\
		{																																								\
			if (TEST::exception == 1) /* dummy to avoid compiler warnings */							\
				TEST::exception++;																													\
			std::cout << std::endl << "    (caught exception of type `"										\
			          << e.getName() << "'";																							\
			if ((e.getLine() > 0) && (std::strcmp(e.getFile(),"")!=0))										\
				std::cout << " outside a subtest, which was thrown in line " << e.getLine()	\
									<< " of file " << e.getFile()																			\
									<< " in function " << e.getFunction();														\
			std::cout << " - unexpected!) " << std::endl;																	\
			std::cout << "    (message is: " << e.what() << ")" << std::endl;							\
		}																																								\
  }																																									\
	/* catch all non-OpenMS exceptions */																							\
	catch (...)																																				\
	{																																									\
		TEST::this_test = false;																												\
		TEST::test = false;																															\
		TEST::all_tests = false;																												\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))						\
		{																																								\
			std::cout << std::endl																												\
			 << "    (caught unidentified and unexpected exception outside a subtest!) "	\
			 << std::endl;																																\
		}																																								\
	}																																									\
	/* check validity of temporary files if known */																	\
	if (!OpenMS::Internal::ClassTest::validate(TEST::tmp_file_list))									\
	{																																									\
		TEST::all_tests = false;																												\
	}																																									\
	/* clean up all temporary files */																								\
	while (TEST::tmp_file_list.size() > 0 && TEST::verbose < 1)												\
	{																																									\
		unlink(TEST::tmp_file_list.back().c_str());																			\
		TEST::tmp_file_list.pop_back();																									\
	}																																									\
	/* check for exit code */																													\
	if (!TEST::all_tests)																															\
	{																																									\
		std::cout << "FAILED";																													\
		if (TEST::add_message != "") std::cout << " (" << TEST::add_message << ")";			\
		std::cout << std::endl;																													\
		return 1;																																				\
	} else {																																					\
		std::cout << "PASSED";																													\
		if (TEST::add_message != "") std::cout << " (" << TEST::add_message << ")";			\
		std::cout << std::endl;																													\
		return 0;																																				\
	}                                                                                 \
}

/**
	@brief Sets an additional text that is displayed after final result of the test.
	
	This can be used to provide additional information about the test to the user.
	It is e.g. used to indicate that the DB test were skipped, when there are no
	credentials given.
	
	@hideinitializer
*/
#define ADD_MESSAGE(message)											\
	TEST::add_message = message;										\

/**	@brief Declare subtest name.

  This macro is used to declare the name of a subtest.  If you want to check
  e.g. the setName method of a class, insert a line #CHECK(setName) in your
  test program. If the test program is called in verbose mode, this leads to
  the name of the subtest being printed on execution.

	This macro also opens a <code>try</code> block to catch any unexpected
	exceptions thrown in the course of a subtest. To catch <em>wanted</em>
	exceptions (i.e. to check for exceptions that are the expected result of
	some command) use the #TEST_EXCEPTION macro.  The <code>try</code> block
	opened by CHECK is closed in #RESULT, so these two macros have to be
	balanced.

	 @hideinitializer
*/
#define CHECK(name_of_test)											\
	TEST::test = true;														\
	TEST::newline = false;												\
	TEST::test_name = #name_of_test;							\
	TEST::test_count = 0;                         \
	TEST::check_line = __LINE__;									\
	if (TEST::verbose > 0)												\
	std::cout << "checking " << TEST::test_name		\
		 << " ... " << std::flush;									\
	try																						\
	{																							\
		while (true)																\
		{

/**
	@brief Print a status message.

	If tests require longer preparations, #STATUS may be used to print some
	intermediated progress messages.  #STATUS uses <code>cout</code> to print
	these messages (in verbose mode only).  The given stream expression
	<code>message</code> is prefixed by the string <code>status:</code> and
	terminated with a newline. All valid operations on a stream may be performed
	in <code>message</code>.

		<b>Example:</b>
		<code>
		STATUS("just calculated x = " << setprecision(10) << x )
		</code>

*/
#define STATUS(message)																					\
					if (TEST::verbose > 1)																\
					{																											\
						if (!TEST::newline)																	\
						{																										\
							TEST::newline = true;															\
							std::cout << std::endl;														\
						}																										\
						std::cout << "  status (line " << __LINE__ << "): "	\
						 << message << std::endl;														\
					}

	/// Shorthand for STATUS("ok")
#define OK STATUS("ok")

/**	@brief Check subtest result.

		Each elementary test macro updates an internal variable (<b>TEST</b>, defined by
		#START_TEST ) that holds the state of the current subtest.

		#RESULT prints whether the subtest has failed or passed in verbose mode
		and updates the internal variables <b>TEST::all_tests</b> that describes
		the state of the whole class test. <b>TEST::all_tests</b> is initialized
		to be <b>true</b>.  If any elementary test fails, <b>TEST::test</b>
		becomes <b>false</b>.  At the time of the next call to #RESULT,
		<b>TEST::all_tests</b> will be set to false, if <b>TEST::test</b> is
		false. One failed elementary test leads therefore to a failed subtest,
		which leads to a failed class test.

		This macro closes the <code>try</code> block opened by #CHECK, so #CHECK and
		#RESULT have to be balanced, or some ugly compile-time errors may occur.
		#RESULT first tries to catch all <code>OpenMS</code> exceptions
		(i.e. exceptions derived from Exception::BaseException). If this fails, it tries to
		catch any exception. After the exception is thrown, the execution will
		continue with the next subtest, the current subtest will be marked as
		failed (as is the whole test program).

		 @hideinitializer
*/
#define RESULT																																											\
			break;																																											\
		}																																																\
  }																																																	\
	/* catch FileNotFound exceptions to print out the file name */																		\
	catch (OpenMS::Exception::FileNotFound& e)																												\
	{																																																	\
		TEST::this_test = false;																																				\
		TEST::test = false;																																							\
		TEST::all_tests = false;																																				\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))														\
		{																																																\
			if (TEST::exception == 1) /* dummy to avoid compiler warnings */															\
				TEST::exception++;																																					\
			std::cout << std::endl << "    (caught exception of type `"																		\
			          << e.getName() << "'";																															\
			if ((e.getLine() > 0) && (std::strcmp(e.getFile(),"")!=0))																		\
				std::cout << " outside a subtest, which was thrown in line " << e.getLine()									\
									<< " of file " << e.getFile()																											\
									<< " in function `" << e.getFunction();																						\
			std::cout << " - unexpected!) " << std::endl;																									\
		}																																																\
  }																																																	\
  catch (::OpenMS::Exception::BaseException& e)																															\
  {																																																	\
    TEST::this_test = false;																																				\
    TEST::test = false;																																							\
    TEST::all_tests = false;																																				\
    if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))														\
    {																																																\
			if (!TEST::newline)																																						\
			{																																															\
				TEST::newline = true;																																				\
				std::cout << std::endl;																																			\
			}																																															\
      std::cout << "    (caught exception of type `"																								\
                << e.getName() << "'";																															\
      if ((e.getLine() > 0) && (std::strcmp(e.getFile(),"")!=0))																		\
				std::cout << " outside a subtest, which was thrown in line " << e.getLine()									\
									<< " of file " << e.getFile()																											\
									<< " in function `" << e.getFunction();																						\
      std::cout << "' - unexpected!) " << std::endl;																								\
			std::cout << "    (message is: `" << e.what() << "')" << std::endl;														\
    }																																																\
  }																																																	\
  catch (std::exception& e)																																					\
  {																																																	\
    TEST::this_test = false;																																				\
    TEST::test = false;																																							\
    TEST::all_tests = false;																																				\
    if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))														\
    {																																																\
			if (!TEST::newline)																																						\
			{																																															\
				TEST::newline = true;																																				\
				std::cout << std::endl;																																			\
			}																																															\
      std::cout << "    (caught std::exception. Cause: `" << e.what() << "')" << std::endl;					\
    }																																																\
  }																																																	\
  catch (std::string& e)																																						\
  {																																																	\
    TEST::this_test = false;																																				\
    TEST::test = false;																																							\
    TEST::all_tests = false;																																				\
    if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))														\
    {																																																\
			if (!TEST::newline)																																						\
			{																																															\
				TEST::newline = true;																																				\
				std::cout << std::endl;																																			\
			}																																															\
      std::cout << "    (caught std::string as an exception: `" << e << "')" << std::endl;					\
    }																																																\
  }																																																	\
  catch (const char* e)																																							\
  {																																																	\
    TEST::this_test = false;																																				\
    TEST::test = false;																																							\
    TEST::all_tests = false;																																				\
    if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))														\
    {																																																\
			if (!TEST::newline)																																						\
			{																																															\
				TEST::newline = true;																																				\
				std::cout << std::endl;																																			\
			}																																															\
      std::cout << "    (caught char pointer as an exception: `" << e << "')" << std::endl;					\
    }																																																\
  }																																																	\
  catch (...)																																												\
  {																																																	\
    TEST::this_test = false;																																				\
    TEST::test = false;																																							\
    TEST::all_tests = false;																																				\
    if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))														\
    {																																																\
			if (!TEST::newline)																																						\
			{																																															\
				TEST::newline = true;																																				\
				std::cout << std::endl;																																			\
			}																																															\
      std::cout << "    (caught unidentified and unexpected exception!)" << std::endl;							\
    }																																																\
  }																																																	\
																																																		\
	TEST::all_tests = TEST::all_tests && TEST::test;																									\
	if (TEST::verbose > 0)																																						\
	{																																																	\
		if (TEST::test){																																								\
			std::cout << (TEST::verbose > 1 ? "passed\n" : "passed") << std::endl;												\
		} else {																																												\
			if (TEST::verbose > 1 )																																				\
			{																																															\
				std::cout << "############################################################\n";							\
			}																																															\
			std::cout																																											\
				<< __FILE__ ":" << TEST::check_line << ":  CHECK(" << TEST::test_name << ")  failed"				\
			  "\n" __FILE__ ":" << __LINE__ << ":  RESULT == failed" << (TEST::verbose > 1 ? "\n" : "" )	\
				<< std::endl;																																								\
		}																																																\
	}                                                                                     \
	/* issue a warning if no tests were performed (unless in destructor)*/ \
	if (TEST::test_count==0)                                                                \
	{                                                                                 \
		bool destructor = false;                                                               \
		for (unsigned int i=0;i!=TEST::test_name.size();++i)                                  \
		{                                                                                     \
			if (TEST::test_name[i] == '~')                                                     \
			{                                                                                 \
				destructor = true;                                                             \
				break;                                                                      \
			}                                                                                \
		}                                                                                \
		if (!destructor) std::cerr << "Warning: no subtests performed in '" << TEST::test_name << "' (line " << __LINE__ << ")!" << std::endl;\
	}																																						      \


/** @brief Create a temporary filename.

  This macro assigns a new temporary filename to the string variable given as
  its argument. The filename is created using the filename of the test and the
  line number where this macro is invoked, for example 'Matrix_test.C' might
  create a temporary file 'Matrix_test_268.tmp' if NEW_TMP_FILE is used in
  line 268.  All temporary files are deleted if #END_TEST is called.  @param
  filename string will contain the filename on completion of the macro.
  
  All temporary files are validated using the XML schema,if the type of file 
  can be determined by FileHandler. Therefor for each file written in a test
  NEW_TMP_FILE should be called. Otherwise only the last writen file is checked.

	@hideinitializer
*/
#define NEW_TMP_FILE(filename)																												\
					filename = OpenMS::Internal::ClassTest::tmpFileName(__FILE__,__LINE__);				\
					TEST::tmp_file_list.push_back(filename);																		\
					if (TEST::verbose > 1)																											\
					{																																						\
						if (!TEST::newline)																												\
						{																																					\
							TEST::newline = true;																										\
							std::cout << std::endl;																									\
						}																																					\
						std::cout << "  creating new temporary filename '" << filename						\
							<< "' (line " << __LINE__ << ")" << std::endl;													\
						std::cout << filename	 << ":0:  [start of new tmp file]" << std::endl;		\
					}



/**	@brief Floating point equality macro.

  Checks whether the absolute value of the difference of the two floating
  point values <b>a</b> and <b>b</b> is less or equal to the value defined by
  #PRECISION.

  @note This macro evaluates its arguments once or twice, depending on verbosity settings.

  @param a floating point value to test
	@param b expected value

   @hideinitializer
*/
#define TEST_REAL_EQUAL(a,b)																																		\
	++TEST::test_count;                                                                           \
	TEST::test_line = __LINE__;																																		\
	TEST::this_test = (fabs((double)(a) -  (double)(b)) < TEST::precision);												\
	TEST::test = TEST::test && TEST::this_test;																										\
	if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))													\
	{																																															\
		if (!TEST::newline)																																					\
		{																																														\
			TEST::newline = true;																																			\
			std::cout << std::endl;																																		\
		}																																														\
		if (TEST::this_test)																																				\
			std::cout << "    (line " << __LINE__ << ":  TEST_REAL_EQUAL(" #a "," #b "): got "				\
				<< (a) << ", expected " << (b) << ")    + " << std::endl;																\
		else																																												\
			std::cout << __FILE__ ":" << TEST::test_line << ":  TEST_REAL_EQUAL(" #a "," #b "): got "	\
				<< (a) << ", expected " << (b) << " (difference is "<< (double)(a) -  (double)(b) << ")    - " << std::endl;																\
	}

/**	@brief String equality macro.

  Both arguments are converted to std::string and tested for equality.
	(That is, we check whether <code>(std::string(a) == std::string(b))</code> holds.)

  @note This macro evaluates its arguments once or twice, depending on verbosity settings.

  @param a value to test
	@param b expected value

   @hideinitializer
*/
#define TEST_STRING_EQUAL(a,b)                                                                           \
	++TEST::test_count;																																				          \
	TEST::test_line = __LINE__;																																					\
	TEST::this_test = (std::string(a) == std::string(b));																								\
	TEST::test = TEST::test && TEST::this_test;																													\
	if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))																\
	{																																																		\
		if (!TEST::newline)																																								\
		{																																																	\
			TEST::newline = true;																																						\
			std::cout << std::endl;																																					\
		}																																																	\
		if (TEST::this_test)																																							\
			std::cout << "    (line " << __LINE__ << " TEST_STRING_EQUAL("<< #a << "," << #b << "): got \""	\
				<< (a) << "\", expected \"" << (b) << "\")    + " << std::endl;																\
		else																																															\
			std::cout << __FILE__ ":" << TEST::test_line << ":  TEST_STRING_EQUAL(" #a "," #b "): got \""		\
				<< (a) << "\", expected \"" << (b) << "\"    - " << std::endl;																\
	}

/**	@brief Generic equality macro.

  This macro uses the operator == to check its two arguments for equality.
  Besides handling some internal stuff, it basically evaluates #((a) == (b))#.

	Remember that operator == has to be defined somehow for the two argument
	types.

  @note This macro evaluates its arguments once or twice, depending on verbosity settings.

	@param a value/object to test
	@param b expected value

	 @hideinitializer
*/
#define TEST_EQUAL(a,b)																																				\
	{                                                                                             \
		++TEST::test_count;																																				\
		TEST::test_line = __LINE__;																																\
		TEST::this_test = ((a) == (b));																														\
		TEST::test = TEST::test && TEST::this_test;																								\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))											\
		{																																													\
			if (!TEST::newline)																																			\
			{																																												\
				TEST::newline = true;																																	\
				std::cout << std::endl;																																\
			}																																												\
			if (TEST::this_test)																																		\
				std::cout << "    (line " << __LINE__ << ":  TEST_EQUAL(" #a "," #b "): got "					\
					<< (a) << ", expected " << (b) << ")    + " << std::endl;														\
			else																																										\
				std::cout << __FILE__ ":" << TEST::test_line << ":  TEST_EQUAL(" #a "," #b "): got "	\
					<< (a) << ", expected " << (b) << "    - " << std::endl;														\
		}																																													\
	}

/**	@brief Generic inequality macro.

  This macro checks for inequality as #TEST_EQUAL tests for equality.  The
  only difference between the two macros is that<b> TEST_NOT_EQUAL</b>
  evaluates #!((a) == (b))#.


	@param a value/object to test
	@param b forbidden value

	 @hideinitializer
*/
#define TEST_NOT_EQUAL(a,b)																																				\
	{                                                                                              \
		++TEST::test_count;																																					\
		TEST::test_line = __LINE__;																																		\
		TEST::this_test = !((a) == (b));																															\
		TEST::test = TEST::test && TEST::this_test;																										\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))													\
		{																																															\
			if (!TEST::newline)																																					\
			{																																														\
				TEST::newline = true;																																			\
				std::cout << std::endl;																																		\
			}																																														\
			if (TEST::this_test)																																				\
			std::cout << "    (line " << __LINE__ << " TEST_NOT_EQUAL(" #a "," #b "): got "							\
				<< (a) << ", forbidden is " << (b) << ")    + " << std::endl;															\
			else																																												\
				std::cout << __FILE__ ":" << TEST::test_line << ":  TEST_NOT_EQUAL(" #a "," #b "): got "	\
				<< (a) << ", forbidden is " << (b) << "    - " << std::endl;															\
		}																																															\
	}

/**	@brief Exception test macro.

  This macro checks if a given type of exception occured while executing the
  given command.  Example: \par #TEST_EXCEPTION(Exception::IndexOverflow,
  vector3[-1])# \par If no or a wrong exception occured, false is returned,
  otherwise true.  @param exception_type the exception-class @param command
  any general C++ or OpenMS-specific command

	 @hideinitializer
*/
#define TEST_EXCEPTION(exception_type,command)																									\
	{                                                                                             \
		++TEST::test_count;																																				\
		TEST::test_line = __LINE__;																																	\
		TEST::exception = 0;																																				\
		try																																													\
		{																																														\
			command;																																									\
		}																																														\
		catch (exception_type)																																			\
		{																																														\
			TEST::exception = 1;																																			\
		}																																														\
		catch (::OpenMS::Exception::BaseException e)																													\
		{																																														\
			TEST::exception = 2;																																			\
			TEST::exception_name = e.getName();																												\
		}																																														\
		catch (...)																																									\
		{																																														\
			TEST::exception = 3;																																			\
		}																																														\
		TEST::this_test = (TEST::exception == 1);																										\
		TEST::test = TEST::test && TEST::this_test;																									\
																																																\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))												\
		{																																														\
			if (!TEST::newline)																																				\
			{																																													\
				TEST::newline = true;																																		\
				std::cout << std::endl;																																	\
			}																																													\
			switch (TEST::exception)																																	\
			{																																													\
				case 0:																																									\
					std::cout << __FILE__ ":" << TEST::test_line <<																				\
					":  TEST_EXCEPTION(" #exception_type "," #command "): no exception thrown!    - "			\
					<< std::endl; break;																																	\
				case 1:																																									\
					std::cout << "    (line " << TEST::test_line <<																				\
					" TEST_EXCEPTION(" #exception_type "," #command "): OK)    +"													\
					<< std::endl; break;																																	\
				case 2:																																									\
					std::cout << __FILE__ ":" << TEST::test_line <<																				\
					":  TEST_EXCEPTION(" #exception_type "," #command "): wrong exception thrown!  \""		\
					<< TEST::exception_name << "\"    - " << std::endl; break;														\
				case 3:																																									\
					std::cout << __FILE__ ":" << TEST::test_line <<																				\
					":  TEST_EXCEPTION(" #exception_type "," #command "): wrong exception thrown!     - "	\
					<< std::endl; break;																																	\
			}																																													\
		}																																														\
	}

/**	
	@brief Macro that suppresses the warning issued when no subtests are performed
	
	Please use this macro only if the method cannot be tested at all or cannot be 
	tested properly on its own. In the later case, the method must however be tested 
	in tests of related methods.
	
	@hideinitializer
*/
#define NOT_TESTABLE                                                                   \
  TEST::test_count = 1;


/**	@brief Exception test macro (with test for exception message).

  This macro checks if a given type of exception occured while executing the
  given command and additionally tests for the message of the exception.
  Example: \par #TEST_EXCEPTION(Exception::IndexOverflow,
  vector3[-1], "a null pointer was specified")# \par If no, a wrong exception occured or a wrong message is
  returned, false is returned, otherwise true.
  @param exception_type the exception-class @param command
  any general C++ or OpenMS-specific command @param message the message
  the exception should give

	 @hideinitializer
*/
#define TEST_EXCEPTION_WITH_MESSAGE(exception_type,command, message)																									\
	{                                                                                            \
		++TEST::test_count;																																				\
		TEST::test_line = __LINE__;																																	\
		TEST::exception = 0;																																				\
		try																																													\
		{																																														\
			command;																																									\
		}																																														\
		catch (exception_type et)																																			\
		{																																														\
			if ( std::string(et.getMessage()) != std::string(message) )																																			\
			{																																			\
				TEST::exception = 4;																																			\
				TEST::exception_message = et.getMessage();																																			\
			}																																			\
			else TEST::exception = 1;																																			\
		}																																														\
		catch (::OpenMS::Exception::BaseException e)																													\
		{																																														\
			TEST::exception = 2;																																			\
			TEST::exception_name = e.getName();																												\
		}																																														\
		catch (...)																																									\
		{																																														\
			TEST::exception = 3;																																			\
		}																																														\
		TEST::this_test = (TEST::exception == 1);																										\
		TEST::test = TEST::test && TEST::this_test;																									\
																																																\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))												\
		{																																														\
			if (!TEST::newline)																																				\
			{																																													\
				TEST::newline = true;																																		\
				std::cout << std::endl;																																	\
			}																																													\
			switch (TEST::exception)																																	\
			{																																													\
				case 0:																																									\
					std::cout << __FILE__ ":" << TEST::test_line <<																				\
					":  TEST_EXCEPTION(" #exception_type "," #command ", " #message "): no exception thrown!    - "			\
					<< std::endl; break;																																	\
				case 1:																																									\
					std::cout << "    (line " << TEST::test_line <<																				\
					" TEST_EXCEPTION(" #exception_type "," #command ", " #message "): OK)    +"													\
					<< std::endl; break;																																	\
				case 2:																																									\
					std::cout << __FILE__ ":" << TEST::test_line <<																				\
					":  TEST_EXCEPTION(" #exception_type "," #command ", " #message "): wrong exception thrown!  \""		\
					<< TEST::exception_name << "\"    - " << std::endl; break;														\
				case 3:																																									\
					std::cout << __FILE__ ":" << TEST::test_line <<																				\
					":  TEST_EXCEPTION(" #exception_type "," #command ", " #message "): wrong exception thrown!     - "	\
					<< std::endl; break;																																	\
				case 4:																																									\
					std::cout << __FILE__ ":" << TEST::test_line <<																				\
					":  TEST_EXCEPTION(" #exception_type "," #command ", " #message "): exception has wrong message: got '"	\
					<< (TEST::exception_message) << "', expected '" << (message) << "'    - "<< std::endl; break;																																	\
			}																																													\
		}																																														\
	}

/** @brief Skip remainder of subtest.

  If the condition is not fulfilled, the remainder of the test is skipped.
  The status (whether it fails or passes) remains unchanged.

   @hideinitializer
*/
#define ABORT_IF(condition)																							\
																																				\
  if (condition)																												\
	{																																			\
		if (TEST::verbose > 1)																							\
		{																																		\
			if (!TEST::newline)																								\
			{																																	\
				TEST::newline = true;																						\
				std::cout << std::endl;																					\
			}																																	\
			std::cout << __FILE__ ":" <<  __LINE__ << ":  ABORT_IF(" #condition "):  TEST ABORTED" \
								<< std::endl;																						\
		}																																		\
		break;																															\
	}

/**
	@brief File comparison macro.

	This macro is used to test file operations. It
	compares the file with name <code>filename</code> against a template file
	<code>templatename</code>. Corresponding lines of the two files have to be
	identical.

	@note line length is limited to 64k characters

  @note This macro evaluates its arguments once or twice, depending on verbosity settings.

	 @hideinitializer
*/
#define TEST_FILE(filename, templatename)																																			\
																																																							\
	{                                                                                                    				\
		++TEST::test_count;																																												\
		TEST::equal_files = true;																																									\
		TEST::infile.open(filename, std::ios::in);																																\
		TEST::templatefile.open(templatename, std::ios::in);																											\
																																																							\
		if (TEST::infile.good() && TEST::templatefile.good())																											\
		{																																																					\
			std::string TEST_FILE__template_line;																																		\
			std::string TEST_FILE__line;																																						\
																																																							\
			while (TEST::infile.good() && TEST::templatefile.good())																								\
			{																																																				\
				TEST::templatefile.getline(TEST::line_buffer, 65535);																									\
				TEST_FILE__template_line = TEST::line_buffer;																													\
				TEST::infile.getline(TEST::line_buffer, 65535);																												\
				TEST_FILE__line = TEST::line_buffer;																																	\
																																																							\
				TEST::equal_files &= (TEST_FILE__template_line == TEST_FILE__line);																		\
				if (TEST_FILE__template_line != TEST_FILE__line)																											\
				{																																																			\
					if (TEST::verbose > 0)																																							\
					{																																																		\
						if (!TEST::newline)																																								\
						{																																																	\
							TEST::newline = true;																																						\
							std::cout << std::endl;																																					\
						}																																																	\
																																																							\
						std::cout << "   TEST_FILE: line mismatch:\n    got:      '"																			\
							<< TEST_FILE__line << "'\n    expected: '" << TEST_FILE__template_line << "'" << std::endl;			\
					}																																																		\
				}																																																			\
			}																																																				\
		} else {																																																	\
			TEST::equal_files = false;																																							\
																																																							\
			if (TEST::verbose > 0)																																									\
			{																																																				\
				if (!TEST::newline)																																										\
				{																																																			\
					TEST::newline = true;																																								\
					std::cout << std::endl;																																							\
				}																																																			\
																																																							\
				std::cout << "    (line " << __LINE__ << ": TEST_FILE(" << #filename << ", " << #templatename ;				\
				std::cout << ") : " << " cannot open file: \"";																												\
				if (!TEST::infile.good())																																							\
				{																																																			\
					std::cout << filename << "\" (input file) ";																												\
				}																																																			\
				if (!TEST::templatefile.good())																																				\
				{																																																			\
					std::cout << templatename << "\" (template file) ";																									\
				}																																																			\
				std::cout << std::endl;																																								\
																																																							\
			}																																																				\
		}																																																					\
		TEST::infile.close();																																											\
		TEST::templatefile.close();																																								\
		TEST::infile.clear();																																											\
		TEST::templatefile.clear();																																								\
																																																							\
		TEST::this_test = TEST::equal_files;																																			\
		TEST::test = TEST::test && TEST::this_test;																																\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))																			\
		{																																																					\
			if (!TEST::newline)																																											\
			{																																																				\
				TEST::newline = true;																																									\
				std::cout << std::endl;																																								\
			}																																																				\
			std::cout << "    (line " << __LINE__ << ": TEST_FILE("<< #filename << ", " << #templatename << "): ";	\
			if (TEST::this_test)																																										\
			{																																																				\
				std::cout << "true + " << std::endl;																																	\
			} else {																																																\
				std::cout << "false - " << std::endl;																																	\
				std::cout << "    (different files:  "<< filename << "  " << templatename << " )\n";									\
			}																																																				\
		}																																																					\
	}

//@} // end of ClassTest

#endif //OPENMS_CONCEPT_CLASSTEST_H
