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

#ifndef OPENMS_CONCEPT_CLASSTEST_H
#define OPENMS_CONCEPT_CLASSTEST_H

//-------------------------------------------------------------------------
//	Note: Emacs has a fantastic command "backslashify" that you can use to
//	line up the backslashes in the define blocks.
//-------------------------------------------------------------------------


/** @brief This indicates that a class test is being compiled.

	Used e.g. in OPENMS_PRECONDITION and OPENMS_POSTCONDITION so that	we
	can	test *these* even if the global OPENMS_DEBUG macro is not set.
*/
#define OPENMS_WITHIN_CLASSTEST 1

// Avoid OpenMS includes here at all costs
// When the included headers are changed, *all* tests have to be recompiled!
// Use the ClassTest class if you need add high-level functionality.
// Includes in the C-file are ok...
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/SYSTEM/ProcessResource.h>

#include <vector>
#include <string>
#include <cstring>
#include <list>
#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>

#ifdef OPENMS_HAS_UNISTD_H
#include <unistd.h> // unlink()
#endif

#include <cstdio>  // tmpnam()
#include <cmath>   // fabs
#include <cstdlib> // getenv()

// Empty declaration to avoid problems in case the namespace is not
// yet defined (e.g. TEST/ClassTest_test.C)

/// Provide a point of redirection for testing the test macros, see ClassTest_test.C
#ifndef std__cout
#define std__cout std::cout
#endif

namespace OpenMS
{
	namespace Internal
	{
		/// Auxilary class for class tests
		namespace OPENMS_DLLAPI ClassTest
		{

			/**
			@brief Validates the given files against the XML schema (if available)
			@return If all files passed the validation
			*/
			bool validate(const std::vector<std::string>& file_names);

			/// Creates a temporary file name from the test name and the line
			std::string tmpFileName(const std::string& file, int line);

			/** @brief Compare floating point numbers using @em absdiff_max_allowed and
			@em ratio_max_allowed.

			Side effects: Updates #fuzzy_message.
			*/
			bool isRealSimilar(double number_1, double number_2);

			/**@brief Compare strings using @em absdiff_max_allowed and @em ratio_max_allowed.

			Side effects: Updates #absdiff, #ratio, #fuzzy_message, #line_num_1_max
			and #line_num_2_max.
			*/
			bool isStringSimilar( const std::string & string_1, const std::string & string_2);
				
			/**@brief Compare files using @em absdiff_max_allowed and @em ratio_max_allowed.

			Side effects: Updates #absdiff, #ratio, #fuzzy_message, #line_num_1_max
			and #line_num_2_max.
			*/
			bool isFileSimilar( const std::string & filename_1, const std::string & filename_2);
				
			/// make sure we have a newline before results from first subtest
			void initialNewline();

			/// print the text, each line gets a prefix, the marked line number gets a special prefix
			void printWithPrefix(const std::string & text, const int marked = -1);

			/// set the whitelist
			void setWhitelist(const char * const /* file */, const int line, const std::string& whitelist);

			/// Maximum ratio of numbers allowed, see #TOLERANCE_RELATIVE.
			extern double ratio_max_allowed;
			
			/// Maximum ratio of numbers observed so far, see #TOLERANCE_RELATIVE.
			extern double ratio_max;
			
			/// Recent ratio of numbers, see #TOLERANCE_RELATIVE.
			extern double ratio;
			
			/// Maximum absolute difference of numbers allowed, see #TOLERANCE_ABSOLUTE.
			extern double absdiff_max_allowed;
			
			/// Maximum difference of numbers observed so far, see #TOLERANCE_ABSOLUTE.
			extern double absdiff_max;
			
			/// Recent absolute difference of numbers, see #TOLERANCE_ABSOLUTE.
			extern double absdiff;

			extern int line_num_1_max;
			extern int line_num_2_max;

			/// Verbosity level ( "-v" is 1 and "-V" is 2 )
			extern int verbose;

			/// Status of the whole test
			extern bool all_tests;

			/// Status of the current subsection
			extern bool test;

			/// Status of last elementary test
			extern bool this_test;

			/// (Used by various macros. Indicates a rough category of the exception being caught.)
			extern int exception;

			/// (Used by various macros.  Stores the "name" of the exception, if applicable.)
			extern std::string exception_name;

			/// (Used by various macros.  Stores the "message" of the exception, if applicable.)
			extern std::string exception_message;

			/// Name of current subsection
			extern std::string test_name;

			/// Line where current subsection started
			extern int start_section_line;

			/// Line of current elementary test
			extern int test_line;

			/// Version string supplied with #START_TEST
			extern const char* version_string;

			/// List of tmp file names (these will be cleaned up, see #NEW_TMP_FILE)
			extern std::vector<std::string> tmp_file_list;

			/// Questionable file tested by #TEST_FILE_EQUAL
			extern std::ifstream infile;

			/// Template (correct) file used by #TEST_FILE_EQUAL
			extern std::ifstream templatefile;

			/// (A variable used by #TEST_FILE_EQUAL)
			extern bool equal_files;

			/// (A buffer for one line from a file. Used by #TEST_FILE_EQUAL.)
			extern char line_buffer[65536];

			/// Counter for the number of elementary tests within the current subsection.
			extern int test_count;

			/// See #ADD_MESSAGE.
			extern std::string add_message;

			/**@brief Last message from a fuzzy comparison.  Written by
			#isRealSimilar(), #isStringSimilar(), #isFileSimilar().  Read by
			#TEST_REAL_SIMILAR, #TEST_STRING_SIMILAR, #TEST_FILE_SIMILAR;
			*/
			extern std::string fuzzy_message;

			/// (Flags whether a new line is in place, depending on context and verbosity setting.  Used by initialNewline() and some macros.)
			extern bool newline;

		}
	}
}

// A namespace alias - apparently these cannot be documented by doxygen (?)
namespace TEST = OpenMS::Internal::ClassTest;

/**
	@defgroup ClassTest Class test macros

	@brief These macros are used by the test programs in the subdirectory
	<code>OpenMS/source/TEST</code>.

	On successful operation the test program will print out the message "PASSED",
	otherwise "FAILED".

	If called with the @b -v option, the test program prints verbose information
	about subsections.

	If called with the @b -V option, the test program prints even more verbose
	information for every elementary test.

  The implementation is done in namespace #OpenMS::Internal::ClassTest.

	@ingroup Concept

*/
//@{


//@name test and subtest
//@{

/**	@brief Begin of the test program for a given class.  @sa #END_TEST.

	The #START_TEST macro defines the start of the test program for a given
	classname.  The classname is printed together with some information when
	calling the test program with any arguments (except for <code>-v</code> or
	<code>-V</code>).

  The second argument version should take the form "$Id:$".  The actual
	version info will then be filled in by Subversion (the revision control
	system).  If it does not, use "svn help propset" to find out how to include
	"Id" in the property "svn:keywords" for the *_test.C file in question.

	The #START_TEST macro should be the first one to call in a test program. It
	opens a global <code>try</code> block to catch any unwanted exceptions.  If
	any of these exceptions occurs, all tests failed.  Exceptions defined by
	%OpenMS (i.e. exception classes derived from Exception::BaseException)
	provide some additional information that is evaluated by the #END_TEST
	macro.  The #END_TEST macro also closes the <code>try</code> block.  This
	<code>try</code> block should never catch an exception!  All exceptions that
	are thrown due to some malfunction in one of the member functions should be
	caught by the <code>try</code> block created by #START_SECTION and
	#END_SECTION .

	 @hideinitializer
*/
#define START_TEST(class_name, version)																					\
int main(int argc, char **argv)																									\
{																																								\
 																																								\
 TEST::version_string = version;																								\
																																								\
	if (argc == 2)																																\
	{																																							\
		if (!strcmp(argv[1], "-v"))																									\
			TEST::verbose = 1;																												\
		if (!strcmp(argv[1], "-V"))																									\
			TEST::verbose = 2;																												\
	};																																						\
																																								\
	if ((argc > 2) || ((argc == 2) && (TEST::verbose == 0)))											\
	{																																							\
		std::cerr																																		\
			<< "This is " << argv[0] << ", the test program for the\n"								\
			<< #class_name " class.\n"																								\
			"\n"																																			\
			"On successful operation it simply returns PASSED,\n"											\
			"otherwise FAILED is printed.\n"																					\
			"If called with an argument of -v,\n"																			\
			"prints detailed information about individual tests.\n"										\
			"Option -V provides verbose information on every subtest." << std::endl;	\
		return 1;																																		\
	}																																							\
																																								\
	if (TEST::verbose > 0)																												\
	{																																							\
		std__cout << "Version: " << TEST::version_string << std::endl;							\
	}																																							\
																																								\
	const char * env_openms_testtimeout;																					\
	env_openms_testtimeout = getenv ("OPENMS_TESTTIMEOUT");												\
	int timeout = 300;																														\
																																								\
	if ( env_openms_testtimeout )																									\
	{																																							\
		try																																					\
		{																																						\
			timeout = boost::lexical_cast<int>(env_openms_testtimeout);								\
		}																																						\
		catch (boost::bad_lexical_cast&)																						\
		{																																						\
			std__cout <<																															\
				"warning: Cannot parse enviroment variable OPENMS_TESTTIMEOUT="					\
				<< env_openms_testtimeout <<																						\
					" (conversion to number failed).  Using timeout: "										\
				<< timeout << ".\n" << std::endl;																				\
		}																																						\
	}																																							\
	OpenMS::ProcessResource::LimitCPUTime(timeout);																\
																																								\
	try {


/**	@brief End of the test program for a class.  @sa #START_TEST.

	The #END_TEST macro implements the correct termination of the test program
	and should therefore be the last macro to call.  It determines the exit code
	based on all previously run subtests and prints out the message "PASSED" or
	"FAILED".  This macro also closes the global <code>try</code> block opened
	by #START_TEST and contains the related <code>catch</code> clauses. If an
	exception is caught here, the test program fails.

	 @hideinitializer
*/
#define END_TEST																												\
	/* global try block */																								\
	}																																			\
	/* catch FileNotFound exceptions to print out the file name */				\
	catch (OpenMS::Exception::FileNotFound e)															\
	{																																			\
		TEST::this_test = false;																						\
		TEST::test = false;																									\
		TEST::all_tests = false;																						\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))	\
		{																																		\
			if (TEST::exception == 1)																					\
				TEST::exception++;																							\
			std__cout << std::endl << "    (caught exception of type `"				\
								<< e.getName() << "'";																	\
			if ((e.getLine() > 0) && (std::strcmp(e.getFile(),"")!=0))				\
				std__cout << " outside a subtest, which was thrown in line " << e.getLine()	\
									<< " of file " << e.getFile()													\
									<< " in function " << e.getFunction();								\
			std__cout << " - unexpected!) " << std::endl;											\
		}																																		\
	}																																			\
	/* catch OpenMS exceptions to retrieve additional information */			\
	catch (OpenMS::Exception::BaseException& e)														\
	{																																			\
		TEST::this_test = false;																						\
		TEST::test = false;																									\
		TEST::all_tests = false;																						\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))	\
		{																																		\
			if (TEST::exception == 1)																					\
				TEST::exception++;																							\
			std__cout << std::endl << "    (caught exception of type `"				\
								<< e.getName() << "'";																	\
			if ((e.getLine() > 0) && (std::strcmp(e.getFile(),"")!=0))				\
				std__cout << " outside a subtest, which was thrown in line " << e.getLine()	\
									<< " of file " << e.getFile()													\
									<< " in function " << e.getFunction();								\
			std__cout << " - unexpected!) " << std::endl;											\
			std__cout << "    (message is: " << e.what() << ")" << std::endl;	\
		}																																		\
	}																																			\
	/* catch all non-OpenMS exceptions */																	\
	catch (...)																														\
	{																																			\
		TEST::this_test = false;																						\
		TEST::test = false;																									\
		TEST::all_tests = false;																						\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))	\
		{																																		\
			std__cout << std::endl <<																					\
				"    (caught unidentified and unexpected exception outside a subtest!) " << \
				std::endl;																											\
		}																																		\
	}																																			\
	/* check validity of temporary files if known */											\
	if (!TEST::validate(TEST::tmp_file_list))															\
	{																																			\
		TEST::all_tests = false;																						\
	}																																			\
	/* clean up all temporary files */																		\
	while (TEST::tmp_file_list.size() > 0 && TEST::verbose < 1)						\
	{																																			\
		unlink(TEST::tmp_file_list.back().c_str());													\
		TEST::tmp_file_list.pop_back();																			\
	}																																			\
	/* check for exit code */																							\
	if (!TEST::all_tests)																									\
	{																																			\
		std__cout << "FAILED";																							\
		if (TEST::add_message != "") std__cout << " (" << TEST::add_message << ")";	\
		std__cout << std::endl;																							\
		return 1;																														\
	}																																			\
	else																																	\
	{																																			\
		std__cout << "PASSED";																							\
		if (TEST::add_message != "") std__cout << " (" << TEST::add_message << ")";	\
		std__cout << std::endl;																							\
		return 0;																														\
	}																																			\
}

/**	@brief Begin of a subtest with a given name.  @sa #END_SECTION.

	The #START_SECTION macro is used to declare the name of a subtest.  Use this
	to examine a member function of the class which was specified in
	#START_TEST.  If you want to check e.g. the memFunc() method of a class,
	insert a line #START_SECTION(memFunc()) in your test program. If the test
	program is called in verbose mode, this leads to the name of the subtest
	being printed on execution.

	This macro also opens a <code>try</code> block to catch any unexpected
	exceptions thrown in the course of a subtest. To catch <em>wanted</em>
	exceptions (i.e. to check for exceptions that are the expected result of
	some command) use the #TEST_EXCEPTION macro.  The <code>try</code> block
	opened by #START_SECTION is closed in #END_SECTION, so these two macros have to be
	balanced.

	 @hideinitializer
*/
#define START_SECTION(name_of_test)																			\
	TEST::test = true;																										\
	TEST::newline = false;																								\
	TEST::test_name = #name_of_test;																			\
	TEST::test_count = 0;																									\
	TEST::start_section_line = __LINE__;																	\
	if (TEST::verbose > 0)																								\
		std__cout << "checking " << TEST::test_name << " ... " << std::flush;	\
	try																																		\
	{																																			\
		while (true)																												\
		{

/**	@brief End of a subtest.  @sa #START_SECTION.

	The #END_SECTION macro defines the end of a subtest.

  Each elementary test macro updates an internal variable (TEST::test) that
	holds the state of the current subtest.  #END_SECTION prints whether the
	subtest has passed or failed (in verbose mode) and updates the internal
	variables <b>TEST::all_tests</b> that describes the state of the whole class
	test. <b>TEST::all_tests</b> is initialized to be <b>true</b>.  If any
	elementary test fails, <b>TEST::test</b> becomes <b>false</b>.  At the time
	of the next call to #END_SECTION, <b>TEST::all_tests</b> will be set to
	false, if <b>TEST::test</b> is false.  One failed elementary test leads
	therefore to a failed subtest, which leads to a failed class test.

	This macro closes the <code>try</code> block opened by #START_SECTION, so
	#START_SECTION and #END_SECTION have to be balanced, or some ugly
	compile-time errors will occur.  #END_SECTION first tries to catch all
	<code>OpenMS::Exception</code>s (i.e. exceptions derived from
	OpenMS::Exception::BaseException).  If this fails, it tries to catch any
	exception.  After the exception is caught, the execution will continue with
	the next subtest, but the current subtest is marked as failed (as is the
	whole test program).

	@hideinitializer
*/
#define END_SECTION																																				\
			break;																																							\
		}																																											\
	}																																												\
	/* catch FileNotFound exceptions to print out the file name */													\
	catch (OpenMS::Exception::FileNotFound& e)																							\
	{																																												\
		TEST::this_test = false;																															\
		TEST::test = false;																																		\
		TEST::all_tests = false;																															\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))									\
		{																																											\
			if (TEST::exception == 1) /* dummy to avoid compiler warnings */										\
				TEST::exception++;																																\
			std__cout << std::endl << "    (caught exception of type `"													\
								<< e.getName() << "'";																										\
			if ((e.getLine() > 0) && (std::strcmp(e.getFile(),"")!=0))													\
				std__cout << " outside a subtest, which was thrown in line " << e.getLine()				\
									<< " of file " << e.getFile()																						\
									<< " in function `" << e.getFunction();																	\
			std__cout << " - unexpected!) " << std::endl;																				\
		}																																											\
	}																																												\
	catch (::OpenMS::Exception::BaseException& e)																						\
	{																																												\
		TEST::this_test = false;																															\
		TEST::test = false;																																		\
		TEST::all_tests = false;																															\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))									\
		{																																											\
			TEST::initialNewline();																															\
			std__cout << "    (caught exception of type `"																			\
								<< e.getName() << "'";																										\
			if ((e.getLine() > 0) && (std::strcmp(e.getFile(),"")!=0))													\
				std__cout << " outside a subtest, which was thrown in line " << e.getLine()				\
									<< " of file " << e.getFile()																						\
									<< " in function `" << e.getFunction();																	\
			std__cout << "' - unexpected!) " << std::endl;																			\
			std__cout << "    (message is: `" << e.what() << "')" << std::endl;									\
		}																																											\
	}																																												\
	catch (std::exception& e)																																\
	{																																												\
		TEST::this_test = false;																															\
		TEST::test = false;																																		\
		TEST::all_tests = false;																															\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))									\
		{																																											\
			TEST::initialNewline();																															\
			std__cout << "    (caught std::exception. Cause: `" << e.what() << "')" << std::endl;	\
		}																																											\
	}																																												\
	catch (std::string& e)																																	\
	{																																												\
		TEST::this_test = false;																															\
		TEST::test = false;																																		\
		TEST::all_tests = false;																															\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))									\
		{																																											\
			TEST::initialNewline();																															\
			std__cout << "    (caught std::string as an exception: `" << e << "')" << std::endl; \
		}																																											\
	}																																												\
	catch (const char* e)																																		\
	{																																												\
		TEST::this_test = false;																															\
		TEST::test = false;																																		\
		TEST::all_tests = false;																															\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))									\
		{																																											\
			TEST::initialNewline();																															\
			std__cout << "    (caught char pointer as an exception: `" << e << "')" << std::endl;	\
		}																																											\
	}																																												\
	catch (...)																																							\
	{																																												\
		TEST::this_test = false;																															\
		TEST::test = false;																																		\
		TEST::all_tests = false;																															\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))									\
		{																																											\
			TEST::initialNewline();																															\
			std__cout << "    (caught unidentified and unexpected exception!)" << std::endl;		\
		}																																											\
	}																																												\
			TEST::all_tests = TEST::all_tests && TEST::test;																		\
	if (TEST::verbose > 0)																																	\
	{																																												\
		if (TEST::test)																																				\
		{																																											\
			std__cout << ": passed" << std::endl;																								\
		}																																											\
		else																																									\
		{																																											\
			if (TEST::verbose > 1 )																															\
			{																																										\
				std__cout << "############################################################\n";		\
			}																																										\
			std__cout <<																																				\
				__FILE__ ":" << TEST::start_section_line <<																				\
				": failed  START_SECTION(" << TEST::test_name << ")\n"														\
				__FILE__ ":" << __LINE__ <<																												\
				": failed  END_SECTION" <<																												\
				(TEST::verbose > 1 ? "\n" : "" ) << std::endl;																		\
		}																																											\
	}																																												\
	/* issue a warning if no tests were performed (unless in destructor)*/									\
	if (TEST::test_count==0)																																\
	{																																												\
		bool destructor = false;																															\
		for (unsigned int i=0;i!=TEST::test_name.size();++i)																	\
		{																																											\
			if (TEST::test_name[i] == '~')																											\
			{																																										\
				destructor = true;																																\
				break;																																						\
			}																																										\
		}																																											\
		if (!destructor) std::cerr << "Warning: no subtests performed in '"										\
															 << TEST::test_name << "' (line " << __LINE__ << ")!" << std::endl;	\
	}																																												\
	if ( TEST::verbose > 1 ) std__cout << std::endl;

//@}


/**	@brief Generic equality macro.

	This macro uses the operator == to check its two arguments for equality.
	Besides handling some internal stuff, it basically evaluates ((a) == (b)).

	Remember that operator == has to be defined somehow for the two argument
	types.

	@note This macro evaluates its arguments once or twice, depending on verbosity settings.

	@param a value/object to test
	@param b expected value

	 @hideinitializer
*/
#define TEST_EQUAL(a,b)																									\
	{																																			\
		++TEST::test_count;																									\
		TEST::test_line = __LINE__;																					\
		TEST::this_test = ((a) == (b));																			\
		TEST::test = TEST::test && TEST::this_test;													\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))	\
		{																																		\
			TEST::initialNewline();																						\
			if (TEST::this_test)																							\
				std__cout << "    (line " << __LINE__ <<												\
					":  TEST_EQUAL(" #a "," #b																		\
					"): got " << (a) <<																						\
					", expected " << (b) <<																				\
					")    + " << std::endl;																				\
			else																															\
				std__cout << __FILE__ ":" << TEST::test_line <<									\
					":  TEST_EQUAL(" #a "," #b																		\
					"): got " << (a) <<																						\
					", expected " << (b) <<																				\
					"    - " << std::endl;																				\
		}																																		\
	}

/**	@brief Generic inequality macro.

	This macro checks for inequality just like #TEST_EQUAL tests for equality.
	The only difference between the two macros is that #TEST_NOT_EQUAL evaluates
	!((a) == (b)).

	@param a value/object to test
	@param b forbidden value

	 @hideinitializer
*/
#define TEST_NOT_EQUAL(a,b)																							\
	{																																			\
		++TEST::test_count;																									\
		TEST::test_line = __LINE__;																					\
		TEST::this_test = !((a) == (b));																		\
		TEST::test = TEST::test && TEST::this_test;													\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))	\
		{																																		\
			TEST::initialNewline();																						\
			if (TEST::this_test)																							\
				std__cout << "    (line " << __LINE__ <<												\
					" TEST_NOT_EQUAL(" #a "," #b																	\
					"): got " << (a) <<																						\
					", forbidden is " << (b) <<																		\
					")    + " << std::endl;																				\
			else																															\
				std__cout << __FILE__ ":" << __LINE__ <<												\
					":  TEST_NOT_EQUAL(" #a "," #b																\
					"): got " << (a) <<																						\
					", forbidden is " << (b) <<																		\
					"    - " << std::endl;																				\
		}																																		\
	}

/**	@brief String equality macro.

	Both arguments are converted to std::string and tested for equality.  (That
	is, we check whether <code>(std::string(a) == std::string(b))</code> holds.)

	@note This macro evaluates its arguments once or twice, depending on verbosity settings.

	@param a value to test
	@param b expected value

	@hideinitializer
*/
#define TEST_STRING_EQUAL(a,b)																					\
	{																																			\
		++TEST::test_count;																									\
		TEST::test_line = __LINE__;																					\
		TEST::this_test = (std::string(a) == std::string(b));								\
		TEST::test = TEST::test && TEST::this_test;													\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))	\
		{																																		\
			TEST::initialNewline();																						\
			if (TEST::this_test)																							\
				std__cout << "    (line " << __LINE__ <<												\
					" TEST_STRING_EQUAL("<< #a << "," << #b <<										\
					"): got \"" << (a) <<																					\
					"\", expected \"" << (b) <<																		\
					"\")    + " << std::endl;																			\
			else																															\
				std__cout << __FILE__ ":" << TEST::test_line <<									\
					":  TEST_STRING_EQUAL(" #a "," #b															\
					"): got \"" << (a) <<																					\
					"\", expected \"" << (b) <<																		\
					"\"    - " << std::endl;																			\
		}																																		\
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
#define TEST_FILE_EQUAL(filename, templatename)																						\
	{																																												\
		++TEST::test_count;																																		\
		TEST::equal_files = true;																															\
		TEST::infile.open(filename, std::ios::in);																						\
		TEST::templatefile.open(templatename, std::ios::in);																	\
																																													\
		if (TEST::infile.good() && TEST::templatefile.good())																	\
		{																																											\
			std::string TEST_FILE__template_line;																								\
			std::string TEST_FILE__line;																												\
																																													\
			while (TEST::infile.good() && TEST::templatefile.good())														\
			{																																										\
				TEST::templatefile.getline(TEST::line_buffer, 65535);															\
				TEST_FILE__template_line = TEST::line_buffer;																			\
				TEST::infile.getline(TEST::line_buffer, 65535);																		\
				TEST_FILE__line = TEST::line_buffer;																							\
																																													\
				TEST::equal_files &= (TEST_FILE__template_line == TEST_FILE__line);								\
				if (TEST_FILE__template_line != TEST_FILE__line)																	\
				{																																									\
					if (TEST::verbose > 0)																													\
					{																																								\
						TEST::initialNewline();																												\
						std__cout << "   TEST_FILE_EQUAL: line mismatch:\n    got:      '"						\
											<< TEST_FILE__line << "'\n    expected: '"													\
											<< TEST_FILE__template_line << "'" << std::endl;										\
					}																																								\
				}																																									\
			}																																										\
		}																																											\
		else																																									\
		{																																											\
			TEST::equal_files = false;																													\
			if (TEST::verbose > 0)																															\
			{																																										\
				TEST::initialNewline();																														\
				std__cout << "    (line " << __LINE__ << ": TEST_FILE_EQUAL(" << #filename << ", " << #templatename ;	\
				std__cout << ") : " << " cannot open file: \"";																		\
				if (!TEST::infile.good())																													\
				{																																									\
					std__cout << filename << "\" (input file) ";																		\
				}																																									\
				if (!TEST::templatefile.good())																										\
				{																																									\
					std__cout << templatename << "\" (template file) ";															\
				}																																									\
				std__cout << std::endl;																														\
																																													\
			}																																										\
		}																																											\
		TEST::infile.close();																																	\
		TEST::templatefile.close();																														\
		TEST::infile.clear();																																	\
		TEST::templatefile.clear();																														\
																																													\
		TEST::this_test = TEST::equal_files;																									\
		TEST::test = TEST::test && TEST::this_test;																						\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))									\
		{																																											\
			TEST::initialNewline();																															\
			std__cout << "    (line " << __LINE__ << ": TEST_FILE_EQUAL("<< #filename << ", " << #templatename << "): "; \
			if (TEST::this_test)																																\
			{																																										\
				std__cout << "true + " << std::endl;																							\
			}																																										\
 			else																																								\
			{																																										\
				std__cout << "false - " << std::endl;																							\
				std__cout << "    (different files:  "<< filename << "  " << templatename << " )\n"; \
			}																																										\
		}																																											\
	}

/**	@brief Floating point similarity macro.

	Checks whether the two numbers are sufficiently close based upon the
	settings of #TOLERANCE_ABSOLUTE and #TOLERANCE_RELATIVE.

	@note This macro evaluates its arguments once or twice, depending on verbosity settings.

  @note Both arguments are converted to @c double.  The actual comparison is done
  by isRealSimilar().

	@param a value to test
	@param b expected value

	@hideinitializer
*/
#define TEST_REAL_SIMILAR(a,b)																					\
	{																																			\
		++TEST::test_count;																									\
		TEST::test_line = __LINE__;																					\
		TEST::this_test = TEST::isRealSimilar( (a), (b) );									\
		TEST::test = TEST::test && TEST::this_test;													\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))	\
		{																																		\
			TEST::initialNewline();																						\
			if (TEST::this_test)																							\
			{																																	\
				std__cout << "    (line " << __LINE__ <<												\
					":  TEST_REAL_SIMILAR(" #a "," #b															\
					"): got " << precisionWrapper(a) <<														\
					", expected " << precisionWrapper(b) <<												\
					")    + " << std::endl;																				\
			}																																	\
			else																															\
			{																																	\
				std__cout << __FILE__ ":" << TEST::test_line <<									\
					":  TEST_REAL_SIMILAR(" #a "," #b															\
					"): got " << precisionWrapper(a) <<														\
					", expected " << precisionWrapper(b) <<												\
					" (absolute: " << TEST::absdiff <<														\
					" [" << TEST::absdiff_max_allowed <<													\
					"], relative: " << TEST::ratio <<															\
					" [" << TEST::ratio_max_allowed <<														\
					"], message: \"" << TEST::fuzzy_message <<										\
					"\" )    - " << std::endl;																		\
			}																																	\
		}																																		\
	}

/**	@brief String similarity macro.

	Compares the two strings using @em FuzzyStringComparator with the settings of
	#TOLERANCE_ABSOLUTE and #TOLERANCE_RELATIVE.

	@note This macro evaluates its arguments once or twice, depending on verbosity settings.

  @note Both arguments are converted to @c std::string.  The actual comparison
  is done by isStringSimilar().

	@param a value to test
	@param b expected value

	@hideinitializer
*/
#define TEST_STRING_SIMILAR(a,b)																				\
	{																																			\
		++TEST::test_count;																									\
		TEST::test_line = __LINE__;																					\
		TEST::this_test = TEST::isStringSimilar( (a), (b) );								\
		TEST::test = TEST::test && TEST::this_test;													\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))	\
		{																																		\
			TEST::initialNewline();																						\
			if (TEST::this_test)																							\
			{																																	\
				std__cout << "    (line " << __LINE__ <<												\
					":  TEST_STRING_SIMILAR(" #a "," #b "):  "										\
					"absolute: " << TEST::absdiff <<															\
					" (" << TEST::absdiff_max_allowed <<													\
					"), relative: " << TEST::ratio <<															\
					" (" << TEST::ratio_max_allowed << ")    +\n";								\
				std__cout << "got:\n";																					\
				TEST::printWithPrefix((a),TEST::line_num_1_max);								\
				std__cout << "expected:\n";																			\
				TEST::printWithPrefix((b),TEST::line_num_2_max);								\
			}																																	\
			else																															\
			{																																	\
				std__cout << __FILE__ ":" << TEST::test_line <<									\
					": TEST_STRING_SIMILAR(" #a "," #b ") ...    -\n"							\
					"got:\n";																											\
				TEST::printWithPrefix((a),TEST::line_num_1_max);								\
				std__cout <<"expected:\n";																			\
				TEST::printWithPrefix((b),TEST::line_num_2_max);								\
				std__cout << "message: \n";																			\
				std__cout << TEST::fuzzy_message;																\
			}																																	\
		}																																		\
	}

/**	@brief File similarity macro.

	Compares the two files using @em FuzzyStringComparator with the settings of
	#TOLERANCE_ABSOLUTE and #TOLERANCE_RELATIVE.

	@note This macro evaluates its arguments once or twice, depending on verbosity settings.

  @note The actual comparison is done by isFileSimilar().

	@param a value to test
	@param b expected value

	@hideinitializer
*/
#define TEST_FILE_SIMILAR(a,b)																					\
	{																																			\
		++TEST::test_count;																									\
		TEST::test_line = __LINE__;																					\
		TEST::this_test = TEST::isFileSimilar( (a), (b) );									\
		TEST::test = TEST::test && TEST::this_test;													\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))	\
		{																																		\
			TEST::initialNewline();																						\
			if (TEST::this_test)																							\
			{																																	\
				std__cout << "    (line " << __LINE__ <<												\
					":  TEST_FILE_SIMILAR(" #a "," #b "):  "											\
					"absolute: " << precisionWrapper(TEST::absdiff) <<						\
					" (" << precisionWrapper(TEST::absdiff_max_allowed) <<				\
					"), relative: " << precisionWrapper(TEST::ratio) <<						\
					" (" << precisionWrapper(TEST::ratio_max_allowed) << ")    +\n"; \
				std__cout << "message: \n";																			\
				std__cout << TEST::fuzzy_message;																\
			}																																	\
			else																															\
			{																																	\
				std__cout << __FILE__ ":" << TEST::test_line <<									\
					": TEST_FILE_SIMILAR(" #a "," #b ") ...    -\n";							\
				std__cout << "message: \n";																			\
				std__cout << TEST::fuzzy_message;																\
			}																																	\
		}																																		\
	}

/**	@brief Define the relative tolerance for floating point comparisons.

	@sa #TEST_REAL_SIMILAR, #TEST_STRING_SIMILAR, #TEST_FILE_SIMILAR

	Several macros consider two numbers sufficiently "close" if <b>the ratio of
	the larger and the smaller</b> is bounded by the value supplied by
	#TOLERANCE_RELATIVE.  The default value is @f$ 1 + 10^{-5} @f$.  It is
	possible to redefine the relative tolerance by calling #TOLERANCE_RELATIVE
	with the new value.

	 @hideinitializer
*/
#define TOLERANCE_RELATIVE(a)																						\
	TEST::ratio_max_allowed = (a);																				\
	if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))	\
	{																																			\
		TEST::initialNewline();																							\
	std__cout << "    (line " << __LINE__ <<															\
		":  TOLERANCE_RELATIVE(" << 	TEST::ratio_max_allowed <<						\
		")   (\""#a"\")" << std::endl;																			\
	}

/**	@brief Define the absolute tolerance for floating point comparisons.

	@sa #TEST_REAL_SIMILAR, #TEST_STRING_SIMILAR, #TEST_FILE_SIMILAR

	Several macros consider two numbers sufficiently "close" if <b>the absolute
	difference</b> is bounded by the value supplied by #TOLERANCE_ABSOLUTE.  The
	default value is @f$ 10^{-5} @f$.  It is possible to redefine the absolute
	tolerance by calling #TOLERANCE_ABSOLUTE with the new value.

	 @hideinitializer
*/
#define TOLERANCE_ABSOLUTE(a)																						\
	TEST::absdiff_max_allowed = (a);																			\
	if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))	\
	{																																			\
		TEST::initialNewline();																							\
		std__cout << "    (line " << __LINE__ <<														\
			":  TOLERANCE_ABSOLUTE(" << 	TEST::absdiff_max_allowed	<<				\
			")   (\""#a"\")" << std::endl;																		\
	}

/** @brief Define the whitelist used by #TEST_STRING_SIMILAR and #TEST_FILE_SIMILAR.

	If both lines contain the same element from this list, they are skipped
	over. (See @em FuzzyStringComparator.)
*/
#define WHITELIST(a)														\
	TEST::setWhitelist(__FILE__,__LINE__,(a));

/**	@brief Exception test macro.

	This macro checks if a given type of exception occured while executing the
	given command.  Example: #TEST_EXCEPTION(Exception::IndexOverflow,
	vector[-1]).  If no or a wrong exception occured, false is returned,
	otherwise true. 

	@param exception_type the exception-class
	@param command any general C++ or OpenMS-specific command

	@hideinitializer
*/
#define TEST_EXCEPTION(exception_type,command)													\
	{																																			\
		++TEST::test_count;																									\
		TEST::test_line = __LINE__;																					\
		TEST::exception = 0;																								\
		try																																	\
		{																																		\
			command;																													\
		}																																		\
		catch (exception_type)																							\
		{																																		\
			TEST::exception = 1;																							\
		}																																		\
		catch (::OpenMS::Exception::BaseException e)												\
		{																																		\
			TEST::exception = 2;																							\
			TEST::exception_name = e.getName();																\
		}																																		\
		catch (...)																													\
		{																																		\
			TEST::exception = 3;																							\
		}																																		\
		TEST::this_test = (TEST::exception == 1);														\
		TEST::test = TEST::test && TEST::this_test;													\
																																				\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))	\
		{																																		\
			TEST::initialNewline();																						\
			switch (TEST::exception)																					\
			{																																	\
			case 0:																														\
				std__cout << __FILE__ ":" << TEST::test_line <<									\
					":  TEST_EXCEPTION(" #exception_type "," #command							\
					"): no exception thrown!    - " << std::endl;									\
				break;																													\
			case 1:																														\
				std__cout << "    (line " << TEST::test_line <<									\
					" TEST_EXCEPTION(" #exception_type "," #command								\
					"): OK)    +" << std::endl;																		\
				break;																													\
			case 2:																														\
				std__cout << __FILE__ ":" << TEST::test_line <<									\
					":  TEST_EXCEPTION(" #exception_type "," #command							\
					"): wrong exception thrown!  \""															\
									<< TEST::exception_name << "\"    - " << std::endl;		\
				break;																													\
			case 3:																														\
				std__cout << __FILE__ ":" << TEST::test_line <<									\
					":  TEST_EXCEPTION(" #exception_type "," #command							\
					"): wrong exception thrown!     - " << std::endl;							\
				break;																													\
			}																																	\
		}																																		\
	}

/**	@brief Exception test macro (with test for exception message).

	This macro checks if a given type of exception occured while executing the
	given command and additionally tests for the message of the exception.

	Example:  #TEST_EXCEPTION(Exception::IndexOverflow, vector[-1], "a null pointer was specified")

	If no, a wrong exception occured or a wrong message is returned, false is
	returned, otherwise true.

	@param exception_type the exception-class
	@param command any general C++ or OpenMS-specific command
	@param message the message the exception should give

	@hideinitializer
*/
#define TEST_EXCEPTION_WITH_MESSAGE(exception_type,command, message)		\
	{																																			\
		++TEST::test_count;																									\
		TEST::test_line = __LINE__;																					\
		TEST::exception = 0;																								\
		try																																	\
		{																																		\
			command;																													\
		}																																		\
		catch (exception_type et)																						\
		{																																		\
			if ( std::string(et.getMessage()) != std::string(message) )				\
			{																																	\
				TEST::exception = 4;																						\
				TEST::exception_message = et.getMessage();											\
			}																																	\
			else TEST::exception = 1;																					\
		}																																		\
		catch (::OpenMS::Exception::BaseException e)												\
		{																																		\
			TEST::exception = 2;																							\
			TEST::exception_name = e.getName();																\
		}																																		\
		catch (...)																													\
		{																																		\
			TEST::exception = 3;																							\
		}																																		\
		TEST::this_test = (TEST::exception == 1);														\
		TEST::test = TEST::test && TEST::this_test;													\
																																				\
		if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))	\
		{																																		\
			TEST::initialNewline();																						\
			switch (TEST::exception)																					\
			{																																	\
			case 0:																														\
					std__cout << __FILE__ ":" << TEST::test_line <<								\
						":  TEST_EXCEPTION(" #exception_type "," #command ", " #message \
						"): no exception thrown!    - " << std::endl;								\
					break;																												\
			case 1:																														\
					std__cout << "    (line " << TEST::test_line <<								\
						" TEST_EXCEPTION(" #exception_type "," #command ", " #message \
						"): OK)    +" << std::endl;																	\
					break;																												\
			case 2:																														\
					std__cout << __FILE__ ":" << TEST::test_line <<								\
						":  TEST_EXCEPTION(" #exception_type "," #command ", " #message \
						"): wrong exception thrown!  \"" <<													\
						TEST::exception_name << "\"    - " << std::endl;						\
					break;																												\
			case 3:																														\
					std__cout << __FILE__ ":" << TEST::test_line <<								\
						":  TEST_EXCEPTION(" #exception_type "," #command ", " #message \
						"): wrong exception thrown!     - " << std::endl;						\
					break;																												\
			case 4:																														\
					std__cout << __FILE__ ":" << TEST::test_line <<								\
						":  TEST_EXCEPTION(" #exception_type "," #command ", " #message \
						"): exception has wrong message: got '" <<									\
						TEST::exception_message <<																	\
						"', expected '" <<																					\
						(message) <<																								\
						"'    - "<< std::endl;																			\
					break;																												\
			}																																	\
		}																																		\
	}

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
#define NEW_TMP_FILE(filename)																					\
	{																																			\
		filename = TEST::tmpFileName(__FILE__,__LINE__);										\
		TEST::tmp_file_list.push_back(filename);														\
		if (TEST::verbose > 1)																							\
		{																																		\
			TEST::initialNewline();																						\
			std__cout << "  creating new temporary filename '" << filename		\
								<< "' (line " << __LINE__ << ")" << std::endl;					\
		}																																		\
	}

/** @brief Skip the remainder of the current subtest.

	If the condition is not fulfilled, the remainder of the current subtest is
	skipped over.  The status (whether it fails or passes) remains unchanged.

	 @hideinitializer
*/
#define ABORT_IF(condition)																							\
	if (condition)																												\
	{																																			\
		if (TEST::verbose > 1)																							\
		{																																		\
			TEST::initialNewline();																						\
			std__cout << __FILE__ ":" <<  __LINE__ <<													\
				":  ABORT_IF(" #condition "):  TEST ABORTED" <<									\
				std::endl;																											\
		}																																		\
		break;																															\
	}

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
		STATUS( "just calculated x = " << precisionWrapper(x) )
		</code>

	@hideinitializer
*/
#define STATUS(message)																			\
	if (TEST::verbose > 1)																		\
	{																													\
		TEST::initialNewline();																	\
		std__cout << __FILE__ ":" <<  __LINE__ << ": status:  "	\
							<< message << std::endl;											\
	}

/**
	@brief Sets an additional text that is displayed after final result of the test.

	This can be used to provide additional information about the test to the user.
	It is e.g. used to indicate that the DB test were skipped, when there are no
	credentials given.

	@hideinitializer
*/
#define ADD_MESSAGE(message)										\
	TEST::add_message = message;


/**
	@brief Macro that suppresses the warning issued when no subtests are performed
	
	Please use this macro only if the method cannot be tested at all or cannot be
	tested properly on its own. In the later case, the method must however be tested
	in tests of related methods.  See also @em test_count.

	@hideinitializer
*/
#define NOT_TESTABLE														\
	TEST::test_count = 1;

//@} // end of ClassTest

#endif //OPENMS_CONCEPT_CLASSTEST_H
