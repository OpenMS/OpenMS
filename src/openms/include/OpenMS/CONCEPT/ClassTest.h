// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl, Chris Bielow, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

// Avoid OpenMS includes here at all costs
// When the included headers are changed, *all* tests have to be recompiled!
// Use the ClassTest class if you need add high-level functionality.
// Includes in ClassTest.cpp are ok...
#include <OpenMS/CONCEPT/PrecisionWrapper.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/CONCEPT/MacrosTest.h>
#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/config.h>

#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <type_traits>
// Empty declaration to avoid problems in case the namespace is not
// yet defined (e.g. TEST/ClassTest_test.cpp)

/// Provide a point of redirection for testing the test macros, see ClassTest_test.cpp
#ifndef stdcout
#define stdcout std::cout
#endif

namespace OpenMS
{
  namespace Internal
  {
    /// Namespace for class tests
    namespace ClassTest
    {

      /**
       @brief Validates the given files against the XML schema (if available)
       @return If all files passed the validation
       */
      bool OPENMS_DLLAPI
      validate(const std::vector<std::string>& file_names);

      /// Creates a temporary file name from the test name and the line with the specified extension
      std::string OPENMS_DLLAPI
      createTmpFileName(const std::string& file, int line, const std::string& extension = "");

      /// This overload returns true; @c float is a floating point type.
      inline bool OPENMS_DLLAPI
      isRealType(float)
      {
        return true;
      }

      /// This overload returns true; @c double is a floating point type.
      inline bool OPENMS_DLLAPI
      isRealType(double)
      {
        return true;
      }

      /// This overload returns true; @c long @c double is a floating point type.
      inline bool OPENMS_DLLAPI
      isRealType(long double)
      {
        return true;
      }

      /// This overload returns true; @c ParamValue will be converted to double by #TEST_REAL_SIMILAR.
      inline bool OPENMS_DLLAPI
      isRealType(const ParamValue&)
      {
          return true;
      }

      /// This overload returns true; @c DataValue will be converted to double by #TEST_REAL_SIMILAR.
      inline bool OPENMS_DLLAPI
      isRealType(const DataValue&)
      {
        return true;
      }

      /// This catch-all template returns false; it will be instantiated for non-floating point types.
      template <typename T>
      inline bool
      isRealType(const T&)
      {
        return false;
      }

      /** @brief Compare floating point numbers using @em absdiff_max_allowed and
       @em ratio_max_allowed.

       Side effects: Updates #fuzzy_message.
       */
      void OPENMS_DLLAPI
      testRealSimilar(const char* file, int line, long double number_1,
                      const char* number_1_stringified,
                      bool number_1_is_realtype, Int number_1_written_digits,
                      long double number_2, const char* number_2_stringified,
                      bool /* number_2_is_realtype */, Int number_2_written_digits);

      /// used by testRealSimilar()
      bool OPENMS_DLLAPI
      isRealSimilar(long double number_1, long double number_2);

      /**@brief Compare strings using @em absdiff_max_allowed and @em ratio_max_allowed.

       This is called by the #TEST_STRING_SIMILAR macro.

       Side effects: Updates #absdiff, #ratio, #fuzzy_message, #line_num_1_max
       and #line_num_2_max.
       */
      void OPENMS_DLLAPI
      testStringSimilar(const char* file, int line,
                        const std::string& string_1,
                        const char* string_1_stringified,
                        const std::string& string_2,
                        const char* string_2_stringified);

      /// used by TEST_STRING_EQUAL
      void OPENMS_DLLAPI
      testStringEqual(const char* file, int line,
                      const std::string& string_1,
                      const char* string_1_stringified,
                      const std::string& string_2,
                      const char* string_2_stringified);

      /**@brief Compare files using @em absdiff_max_allowed and @em ratio_max_allowed.

       Side effects: Updates #absdiff, #ratio, #fuzzy_message, #line_num_1_max
       and #line_num_2_max.
       */
      bool OPENMS_DLLAPI
      isFileSimilar(const std::string& filename_1,
                    const std::string& filename_2);

      /// make sure we have a newline before results from first subtest
      void OPENMS_DLLAPI
      initialNewline();

      /// print the text, each line gets a prefix, the marked line number gets a special prefix
      void OPENMS_DLLAPI
      printWithPrefix(const std::string& text, const int marked = -1);

      /**
         @brief Set up some classtest variables as obtained from the 'START_TEST' macro 
                and check that no additional arguments were passed to the test executable.

         @param version A version string, obtained from 'START_TEST(FuzzyStringComparator, "<VERSION>")'
         @param class_name The class under test (used for error messages etc), obtained from 'START_TEST(FuzzyStringComparator, "<VERSION>")'
         @param argc The number of arguments to the main() function of the class test (must be 1; test will quit otherwise)
         @param argv0 Name of the executable (for debug output)
      */
      void OPENMS_DLLAPI mainInit(const char* version, const char* class_name, int argc, const char* argv0);

      /**
        @brief Test if two files are exactly equal (used in TEST_FILE_EQUAL macro)
        
        @param line The line where the macro was called (for reporting)
        @param filename The temp file
        @param templatename The ground truth file
        @param filename_stringified The expression used as the first macro argument
        @param templatename_stringified The expression used as the second macro argument
      */
      void OPENMS_DLLAPI filesEqual(int line, const char* filename, const char* templatename, const char* filename_stringified, const char* templatename_stringified);

      /// removed all temporary files created with the NEW_TMP_FILE macro
      void OPENMS_DLLAPI removeTempFiles();
      
      /// set the whitelist_
      void OPENMS_DLLAPI
      setWhitelist(const char* const /* file */, const int line,
                   const std::string& whitelist);

      /// Maximum ratio of numbers allowed, see #TOLERANCE_RELATIVE.
      extern OPENMS_DLLAPI double ratio_max_allowed;

      /// Maximum ratio of numbers observed so far, see #TOLERANCE_RELATIVE.
      extern OPENMS_DLLAPI double ratio_max;

      /// Recent ratio of numbers, see #TOLERANCE_RELATIVE.
      extern OPENMS_DLLAPI double ratio;

      /// Maximum absolute difference of numbers allowed, see #TOLERANCE_ABSOLUTE.
      extern OPENMS_DLLAPI double absdiff_max_allowed;

      /// Maximum difference of numbers observed so far, see #TOLERANCE_ABSOLUTE.
      extern OPENMS_DLLAPI double absdiff_max;

      /// Recent absolute difference of numbers, see #TOLERANCE_ABSOLUTE.
      extern OPENMS_DLLAPI double absdiff;

      extern OPENMS_DLLAPI int line_num_1_max;
      extern OPENMS_DLLAPI int line_num_2_max;

      /// Verbosity level ( "-v" is 1 and "-V" is 2 )
      extern OPENMS_DLLAPI int verbose;

      /// Status of the whole test
      extern OPENMS_DLLAPI bool all_tests;

      /// Status of the current subsection
      extern OPENMS_DLLAPI bool test;

      /// Status of last elementary test
      extern OPENMS_DLLAPI bool this_test;

      /// (Used by various macros. Indicates a rough category of the exception being caught.)
      extern OPENMS_DLLAPI int exception;

      /// (Used by various macros.  Stores the "name" of the exception, if applicable.)
      extern OPENMS_DLLAPI std::string exception_name;

      /// (Used by various macros.  Stores the "message" of the exception, if applicable.)
      extern OPENMS_DLLAPI std::string exception_message;

      /// Name of current subsection
      extern OPENMS_DLLAPI std::string test_name;

      /// Line where current subsection started
      extern OPENMS_DLLAPI int start_section_line;

      /// Line of current elementary test
      extern OPENMS_DLLAPI int test_line;

      /// Version string supplied with #START_TEST
      extern OPENMS_DLLAPI const char* version_string;

      /// List of tmp file names (these will be cleaned up, see #NEW_TMP_FILE)
      extern OPENMS_DLLAPI std::vector<std::string> tmp_file_list;

      /// List of all failed lines for summary at the end of the test
      extern OPENMS_DLLAPI std::vector<UInt> failed_lines_list;

      /// Questionable file tested by #TEST_FILE_EQUAL
      extern OPENMS_DLLAPI std::ifstream infile;

      /// Template (correct) file used by #TEST_FILE_EQUAL
      extern OPENMS_DLLAPI std::ifstream templatefile;

      /// (A variable used by #TEST_FILE_EQUAL)
      extern OPENMS_DLLAPI bool equal_files;

      /// (A buffer for one line from a file. Used by #TEST_FILE_EQUAL.)
      extern OPENMS_DLLAPI char line_buffer[65536];

      /// Counter for the number of elementary tests within the current subsection.
      extern OPENMS_DLLAPI int test_count;

      /// See #ADD_MESSAGE.
      extern OPENMS_DLLAPI std::string add_message;

      /**@brief Last message from a fuzzy comparison.  Written by
       #isRealSimilar(), #testStringSimilar(), #isFileSimilar().  Read by
       #TEST_REAL_SIMILAR, #TEST_STRING_SIMILAR, #TEST_FILE_SIMILAR;
       */
      extern OPENMS_DLLAPI std::string fuzzy_message;

      /// (Flags whether a new line is in place, depending on context and verbosity setting.  Used by initialNewline() and some macros.)
      extern OPENMS_DLLAPI bool newline;

      template <typename T1, typename T2>
      void
      testEqual(const char* /*file*/, int line, const T1& expression_1,
                const char* expression_1_stringified,
                const T2& expression_2,
                const char* expression_2_stringified)
      {
        ++test_count;
        test_line = line;
        this_test = bool(expression_1 == T1(expression_2)) ;
        test &= this_test;
        {
          initialNewline();
          if (!this_test || verbose > 1)
          {
            stdcout << ' ' << (this_test ? '+' : '-') << "  line " << line << " : TEST_EQUAL(" << expression_1_stringified << ','
                    << expression_2_stringified << "): got '";
            if constexpr (std::is_enum_v<T1> && std::is_enum_v<T2>)
            {
              stdcout << static_cast<int>(expression_1) << "', expected '" << static_cast<int>(expression_2) << "'\n";
            }
            else
            { 
              stdcout << expression_1 << "', expected '" << expression_2 << "'\n";
            }
          }
          if (!this_test)
          {
            failed_lines_list.push_back(line);
          }
        }
      }

      void testTrue(const char* /*file*/, int line, const bool expression_1, const char* expression_1_stringified)
      {
        ++test_count;
        test_line = line;
        this_test = expression_1;
        test &= this_test;
        {
          initialNewline();
          if (this_test)
          {
            if (verbose > 1)
            {
              stdcout << " +  line " << line << ":  TEST_TRUE(" << expression_1_stringified << "): ok\n";
            }
          }
          else
          {
            stdcout << " -  line " << line << ":  TEST_TRUE(" << expression_1_stringified << "): failed\n";
            failed_lines_list.push_back(line);
          }
        }
      }

      void testFalse(const char* /*file*/, int line, const bool expression_1, const char* expression_1_stringified)
      {
        ++test_count;
        test_line = line;
        this_test = !expression_1;
        test &= this_test;
        {
          initialNewline();
          if (this_test)
          {
            if (verbose > 1)
            {
              stdcout << " +  line " << line << ":  TEST_FALSE(" << expression_1_stringified << "): ok\n";
            }
          }
          else
          {
            stdcout << " -  line " << line << ":  TEST_FALSE(" << expression_1_stringified << "): failed\n";
            failed_lines_list.push_back(line);
          }
        }
      }

      template <typename T1, typename T2>
      void
      testNotEqual(const char* /*file*/, int line, const T1& expression_1,
                   const char* expression_1_stringified,
                   const T2& expression_2,
                   const char* expression_2_stringified)
      {
        ++test_count;
        test_line = line;
        this_test = !(expression_1 == T1(expression_2));
        test &= this_test;
        {
          initialNewline();
          if (!this_test || verbose > 1)
          {
            stdcout << ' ' << (this_test ? '+' : '-') << "  line " << line << " : TEST_NOT_EQUAL(" << expression_1_stringified << ','
                    << expression_2_stringified << "): got '";
            if constexpr (std::is_enum_v<T1> && std::is_enum_v<T2>)
            {
              stdcout << static_cast<int>(expression_1) << "', forbidden is '" << static_cast<int>(expression_2) << "'\n";
            }
            else { stdcout << expression_1 << "', expected '" << expression_2 << "'\n"; }
          }
          if (!this_test)
          {
            failed_lines_list.push_back(line);
          }
        }
      }

      
      void OPENMS_DLLAPI printLastException(std::ostream& out);
      
      int OPENMS_DLLAPI endTestPostProcess(std::ostream& out);

      void OPENMS_DLLAPI endSectionPostProcess(std::ostream& out, const int line);
    }
  }
}

// A namespace alias - apparently these cannot be documented using doxygen (?)
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

 To create a test you can use the 'create_test.php' script in %OpenMS/tools/
 (other useful scripts in the same directory - have a look).

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

 The second argument version should take the form "$Id:$" but is currently
 deprecated.  Originally, the SVN revision was annotated by the revision 
 control system.

 The #START_TEST macro should be the first one to call in a test program. It
 opens a global <code>try</code> block to catch any unwanted exceptions.  If
 any of these exceptions occurs, all tests fail.  Exceptions defined by
 %OpenMS (i.e. exception classes derived from Exception::BaseException)
 provide some additional information that is evaluated by the #END_TEST
 macro.  The #END_TEST macro also closes the <code>try</code> block.  This
 <code>try</code> block should never catch an exception!  All exceptions that
 are thrown due to some malfunction in one of the member functions should be
 caught by the <code>try</code> block created by #START_SECTION and
 #END_SECTION .

 @hideinitializer
 */
#define START_TEST(class_name, version)                                                   \
  int main(int argc, char** argv)                                                         \
  {                                                                                       \
    TEST::mainInit(version, #class_name, argc, argv[0]);                                  \
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
#define END_TEST                                                                          \
  /* global try block */                                                                  \
    }                                                                                     \
    catch (...)                                                                           \
    {                                                                                     \
      TEST::printLastException(stdcout);                                                  \
    }                                                                                     \
    return TEST::endTestPostProcess(stdcout);                                             \
  }

/**	@brief Begin of a subtest with a given name.  @sa #END_SECTION.

 The #START_SECTION macro is used to declare the name of a subtest.  Use this
 to examine a member function of the class which was specified in
 #START_TEST.  If you want to check e.g. the memFunc() method of a class,
 insert a line #START_SECTION(memFunc()) in your test program. If the test
 program is called in verbose mode, this leads to the name of the subtest
 being printed on execution.
 If you are testing a non-public method you can use the [EXTRA] statement,
 e.g. #START_SECTION([EXTRA]memFunc())
 to indicate this. Otherwise you will trigger a warning by %OpenMS/tools/checker.php
 due to this unexpected subtest.

 This macro also opens a <code>try</code> block to catch any unexpected
 exceptions thrown in the course of a subtest. To catch <em>wanted</em>
 exceptions (i.e. to check for exceptions that are the expected result of
 some command) use the #TEST_EXCEPTION macro.  The <code>try</code> block
 opened by #START_SECTION is closed in #END_SECTION, so these two macros have to be
 balanced.

 @hideinitializer
 */
#define START_SECTION(name_of_test)                                                       \
  TEST::test = true;                                                                      \
  TEST::newline = false;                                                                  \
  TEST::test_name = # name_of_test;                                                       \
  TEST::test_count = 0;                                                                   \
  TEST::start_section_line = __LINE__;                                                    \
  stdcout << "checking " << TEST::test_name << " ... " << std::flush;                   \
  try                                                                                     \
  {                                                                                       \
    while (true)                                                                          \
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
#define END_SECTION                                                                       \
  break;                                                                                  \
  }                                                                                       \
  }                                                                                       \
  catch (...)  \
  {   \
    TEST::printLastException(stdcout);\
  } \
  TEST::endSectionPostProcess(stdcout, __LINE__);


//@}

/**	@brief Generic equality macro.

 This macro uses the operator == to check its two arguments for equality.
 Besides handling some internal stuff, it basically evaluates ((a) == (b)).

 Remember that operator == has to be defined somehow for the two argument
 types. Additionally the << operator needs to be defined.
 If only == is available you will get a compilation error. As workaround
 use TEST_EQUAL(a==b, true) thereby making bug tracing harder,
 as you won't see the values of a and b.


 @note This macro evaluates its arguments once or twice, depending on verbosity settings.

 @param a value/object to test
 @param b expected value

 @hideinitializer
 */
#define TEST_EQUAL(a, b) TEST::testEqual(__FILE__, __LINE__, (a), (# a), (b), (# b));

/**	@brief Boolean test macro.

 This macro tests if its argument evaluates to 'true'.
 If possible use TEST_EQUAL(a, b) instead of TEST_TRUE(a==b), because the latter makes bug tracing harder.

 @param a value/object convertible to bool
 
 @hideinitializer
*/
#define TEST_TRUE(a) TEST::testTrue(__FILE__, __LINE__, (a), (#a));

/**	@brief Boolean test macro.

 This macro tests if its argument evaluates to 'false'.
 If possible use TEST_NOT_EQUAL(a, b) instead of TEST_FALSE(a!=b), because the latter makes bug tracing harder.

 @param a value/object convertible to bool

 @hideinitializer
*/
#define TEST_FALSE(a) TEST::testFalse(__FILE__, __LINE__, (a), (#a));


/**	@brief Generic inequality macro.

 This macro checks for inequality just like #TEST_EQUAL tests for equality.
 The only difference between the two macros is that #TEST_NOT_EQUAL evaluates
 !((a) == (b)).

 @param a value/object to test
 @param b forbidden value

 @hideinitializer
 */
#define TEST_NOT_EQUAL(a, b) TEST::testNotEqual(__FILE__, __LINE__, (a), (# a), (b), (# b));

/**	@brief String equality macro.

 Both arguments are converted to std::string and tested for equality.  (That
 is, we check whether <code>(std::string(a) == std::string(b))</code> holds.)

 @note This macro evaluates its arguments once or twice, depending on verbosity settings.

 @param a value to test
 @param b expected value

 @hideinitializer
 */
#define TEST_STRING_EQUAL(a, b) TEST::testStringEqual(__FILE__, __LINE__, (a), (# a), (b), (# b));

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
#define TEST_FILE_EQUAL(filename, templatename)                                           \
  {                                                                                       \
    TEST::filesEqual(__LINE__, filename, templatename, #filename, #templatename);         \
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
#define TEST_REAL_SIMILAR(a, b) TEST::testRealSimilar(__FILE__, __LINE__, (a), (# a), TEST::isRealType(a), writtenDigits(a), (b), (# b), TEST::isRealType(b), writtenDigits(b));

/**	@brief String similarity macro.

 Compares the two strings using @em FuzzyStringComparator with the settings of
 #TOLERANCE_ABSOLUTE and #TOLERANCE_RELATIVE.

 @note This macro evaluates its arguments once or twice, depending on verbosity settings.

 @note Both arguments are converted to @c std::string.  The actual comparison
 is done by testStringSimilar().

 @param a value to test
 @param b expected value

 @hideinitializer
 */
#define TEST_STRING_SIMILAR(a, b) TEST::testStringSimilar(__FILE__, __LINE__, (a), (# a), (b), (# b));

/**	@brief File similarity macro.

 Compares the two files using @em FuzzyStringComparator with the settings of
 #TOLERANCE_ABSOLUTE and #TOLERANCE_RELATIVE.

 @note This macro evaluates its arguments once or twice, depending on verbosity settings.

 @note The actual comparison is done by isFileSimilar().

 @param a value to test
 @param b expected value

 @hideinitializer
 */
#define TEST_FILE_SIMILAR(a, b)                                                           \
  {                                                                                       \
    ++TEST::test_count;                                                                   \
    TEST::test_line = __LINE__;                                                           \
    TEST::this_test = TEST::isFileSimilar((a), (b));                                      \
    TEST::test = TEST::test && TEST::this_test;                                           \
    {                                                                                     \
      TEST::initialNewline();                                                             \
      if (TEST::this_test)                                                                \
      {                                                                                   \
        if (TEST::verbose > 1)                                                            \
        {                                                                                 \
          stdcout << " +  line " << __LINE__                                              \
                    << ":  TEST_FILE_SIMILAR(" # a "," # b "):  absolute: "               \
                    << precisionWrapper(TEST::absdiff)                                    \
                    << " ("                                                               \
                    << precisionWrapper(TEST::absdiff_max_allowed)                        \
                    << "), relative: "                                                    \
                    << precisionWrapper(TEST::ratio)                                      \
                    << " ("                                                               \
                    << precisionWrapper(TEST::ratio_max_allowed)                          \
                    << ")\n";                                                             \
          stdcout << "message: \n";                                                       \
          stdcout << TEST::fuzzy_message;                                                 \
        }                                                                                 \
      }                                                                                   \
      else                                                                                \
      {                                                                                   \
        stdcout << " -  line " << TEST::test_line <<                                      \
          ": TEST_FILE_SIMILAR(" # a "," # b ") ...    -\n";                              \
        stdcout << "message: \n";                                                         \
        stdcout << TEST::fuzzy_message;                                                   \
        TEST::failed_lines_list.push_back(TEST::test_line);                               \
      }                                                                                   \
    }                                                                                     \
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
#define TOLERANCE_RELATIVE(a)                                                             \
  TEST::ratio_max_allowed = (a);                                                          \
  {                                                                                       \
    TEST::initialNewline();                                                               \
    if (TEST::verbose > 1)                                                                \
    {                                                                                     \
      stdcout << " +  line " << __LINE__ <<                                               \
        ":  TOLERANCE_RELATIVE(" <<     TEST::ratio_max_allowed <<                        \
        ")   (\"" # a "\")\n";                                                            \
    }                                                                                     \
  }

/**	@brief Define the absolute tolerance for floating point comparisons.

 @sa #TEST_REAL_SIMILAR, #TEST_STRING_SIMILAR, #TEST_FILE_SIMILAR

 Several macros consider two numbers sufficiently "close" if <b>the absolute
 difference</b> is bounded by the value supplied by #TOLERANCE_ABSOLUTE.  The
 default value is @f$ 10^{-5} @f$.  It is possible to redefine the absolute
 tolerance by calling #TOLERANCE_ABSOLUTE with the new value.

 @hideinitializer
 */
#define TOLERANCE_ABSOLUTE(a)                                                             \
  TEST::absdiff_max_allowed = (a);                                                        \
  {                                                                                       \
    TEST::initialNewline();                                                               \
    if (TEST::verbose > 1)                                                                \
    {                                                                                     \
      stdcout << " +  line " << __LINE__ <<                                               \
        ":  TOLERANCE_ABSOLUTE(" <<     TEST::absdiff_max_allowed   <<                    \
        ")   (\"" # a "\")\n";                                                            \
    }                                                                                     \
  }

/** @brief Define the whitelist_ used by #TEST_STRING_SIMILAR and #TEST_FILE_SIMILAR.

 If both lines contain the same element from this list, they are skipped
 over. (See @em FuzzyStringComparator.)
 */
#define WHITELIST(a) TEST::setWhitelist(__FILE__, __LINE__, (a));

/**	@brief Exception test macro.

 This macro checks if a given type of exception occurred while executing the
 given command.  Example: #TEST_EXCEPTION(Exception::IndexOverflow,
 vector[-1]).  If no or a wrong exception occurred, false is returned,
 otherwise true.

 @param exception_type the exception-class
 @param command any general C++ or OpenMS-specific command

 @hideinitializer
 */
#define TEST_EXCEPTION(exception_type, command)                                           \
  {                                                                                       \
    ++TEST::test_count;                                                                   \
    TEST::test_line = __LINE__;                                                           \
    TEST::exception = 0;                                                                  \
    try                                                                                   \
    {                                                                                     \
      command;                                                                            \
    }                                                                                     \
    catch (exception_type&)                                                               \
    {                                                                                     \
      TEST::exception = 1;                                                                \
    }                                                                                     \
    catch (::OpenMS::Exception::BaseException& e)                                         \
    {                                                                                     \
      TEST::exception = 2;                                                                \
      TEST::exception_name = e.getName();                                                 \
    }                                                                                     \
    catch (const std::exception& e)                                                       \
    {                                                                                     \
      TEST::exception = 3;                                                                \
      TEST::exception_name = e.what();                                                    \
    }                                                                                     \
    catch (...)                                                                           \
    {                                                                                     \
      TEST::exception = 4;                                                                \
    }                                                                                     \
    TEST::this_test = (TEST::exception == 1);                                             \
    TEST::test = TEST::test && TEST::this_test;                                           \
    {                                                                                     \
      TEST::initialNewline();                                                             \
      switch (TEST::exception)                                                            \
      {                                                                                   \
      case 0:                                                                             \
        stdcout << " -  line " << TEST::test_line <<                                    \
          ":  TEST_EXCEPTION(" # exception_type "," # command                               \
          "): no exception thrown!\n";                                                    \
        TEST::failed_lines_list.push_back(TEST::test_line);                               \
        break;                                                                            \
      case 1:                                                                             \
        if (TEST::verbose > 1)                                                            \
        {                                                                                 \
          stdcout << " +  line " << TEST::test_line <<                                    \
            ":  TEST_EXCEPTION(" # exception_type "," # command                           \
            "): OK\n";                                                                    \
        }                                                                                 \
        break;                                                                            \
      case 2:                                                                             \
        stdcout << " -  line " << TEST::test_line <<                                    \
          ":  TEST_EXCEPTION(" # exception_type "," # command                               \
          "): wrong exception thrown!  \""                                                  \
                  << TEST::exception_name << "\"\n";                                     \
        TEST::failed_lines_list.push_back(TEST::test_line);                               \
        break;                                                                            \
      case 3:                                                                             \
        stdcout << " -  line " << TEST::test_line <<                                    \
          ":  TEST_EXCEPTION(" # exception_type "," # command                               \
          "): wrong exception thrown!  \""                                                  \
                  << TEST::exception_name << "\"\n";                                     \
        TEST::failed_lines_list.push_back(TEST::test_line);                               \
        break;                                                                            \
      case 4:                                                                             \
        stdcout << " -  line " << TEST::test_line <<                                    \
          ":  TEST_EXCEPTION(" # exception_type "," # command                               \
          "): wrong exception thrown!\n";                                                 \
        TEST::failed_lines_list.push_back(TEST::test_line);                               \
        break;                                                                            \
      }                                                                                   \
    }                                                                                     \
  }
  
/** @brief Precondition test macro

  This macro checks if a precondition violation is detected while executing the command,
  similar to <code>TEST_EXCEPTION(Exception::Precondition,command)</code>.
  However the test is executed only when the #OPENMS_PRECONDITION macros are active,
  i.e., when compiling in Debug mode.  (See Macros.h)

 @param command any general C++ or OpenMS-specific command

  @hideinitializer
 */
#ifdef OPENMS_ASSERTIONS
#define TEST_PRECONDITION_VIOLATED(command) TEST_EXCEPTION(Exception::Precondition, command);
#else
#define TEST_PRECONDITION_VIOLATED(command) STATUS("TEST_PRECONDITION_VIOLATED(" # command ")  -  skipped");
#endif

/** @brief Postcondition test macro

  This macro checks if a postcondition violation is detected while executing the command,
  similar to <code>TEST_EXCEPTION(Exception::Postcondition,command)</code>.
  However the test is executed only when the #OPENMS_POSTCONDITION macros are active,
  i.e., when compiling in Debug mode.  (See Macros.h)

 @param command any general C++ or OpenMS-specific command

  @hideinitializer
 */
#ifdef OPENMS_ASSERTIONS
#define TEST_POSTCONDITION_VIOLATED(command) TEST_EXCEPTION(Exception::Postcondition, command);
#else
#define TEST_POSTCONDITION_VIOLATED(command) STATUS("TEST_POSTCONDITION_VIOLATED(" # command ")  -  skipped");
#endif


/**	@brief Exception test macro (with test for exception message).

 This macro checks if a given type of exception occurred while executing the
 given command and additionally tests for the message of the exception.

 Example:  #TEST_EXCEPTION_WITH_MESSAGE(Exception::IndexOverflow, vector[-1], "a null pointer was specified")

 If no, a wrong exception occurred or a wrong message is returned, false is
 returned, otherwise true.

 @param exception_type the exception-class
 @param command any general C++ or OpenMS-specific command
 @param message the message the exception should give

 @hideinitializer
 */
#define TEST_EXCEPTION_WITH_MESSAGE(exception_type, command, message)                     \
  {                                                                                       \
    ++TEST::test_count;                                                                   \
    TEST::test_line = __LINE__;                                                           \
    TEST::exception = 0;                                                                  \
    try                                                                                   \
    {                                                                                     \
      command;                                                                            \
    }                                                                                     \
    catch (exception_type& et)                                                            \
    {                                                                                     \
      if (std::string(et.what()) != std::string(message))                                 \
      {                                                                                   \
        TEST::exception = 4;                                                              \
        TEST::exception_message = et.what();                                              \
      }                                                                                   \
      else TEST::exception = 1;                                                           \
    }                                                                                     \
    catch (::OpenMS::Exception::BaseException& e)                                         \
    {                                                                                     \
      TEST::exception = 2;                                                                \
      TEST::exception_name = e.getName();                                                 \
    }                                                                                     \
    catch (...)                                                                           \
    {                                                                                     \
      TEST::exception = 3;                                                                \
    }                                                                                     \
    TEST::this_test = (TEST::exception == 1);                                             \
    TEST::test = TEST::test && TEST::this_test;                                           \
                                                                                          \
    {                                                                                     \
      TEST::initialNewline();                                                             \
      switch (TEST::exception)                                                            \
      {                                                                                   \
      case 0:                                                                             \
        stdcout << " -  line " << TEST::test_line <<                                    \
          ":  TEST_EXCEPTION_WITH_MESSAGE(" # exception_type "," # command ", " # message   \
          "): no exception thrown!\n";                                                    \
        TEST::failed_lines_list.push_back(TEST::test_line);                               \
        break;                                                                            \
      case 1:                                                                             \
        if (TEST::verbose > 1)                                                            \
        {                                                                                 \
          /* this is actually what we want to get:  */                                      \
          stdcout << " +  line " << TEST::test_line <<                                    \
            ":  TEST_EXCEPTION_WITH_MESSAGE(" # exception_type "," # command ", " # message   \
            "): OK\n";                                                                      \
        }                                                                                 \
        break;                                                                            \
      case 2:                                                                             \
        stdcout << " -  line " << TEST::test_line <<                                    \
          ":  TEST_EXCEPTION_WITH_MESSAGE(" # exception_type "," # command ", " # message   \
          "): wrong exception thrown!  \"" <<                                               \
          TEST::exception_name << "\"\n";                                                \
        TEST::failed_lines_list.push_back(TEST::test_line);                               \
        break;                                                                            \
      case 3:                                                                             \
        stdcout << " -  line " << TEST::test_line <<                                    \
          ":  TEST_EXCEPTION_WITH_MESSAGE(" # exception_type "," # command ", " # message   \
          "): wrong exception thrown!\n";                                                 \
        TEST::failed_lines_list.push_back(TEST::test_line);                               \
        break;                                                                            \
      case 4:                                                                             \
        stdcout << " -  line " << TEST::test_line <<                                    \
          ":  TEST_EXCEPTION_WITH_MESSAGE(" # exception_type "," # command ", " # message   \
          "): exception has wrong message: got '" <<                                        \
          TEST::exception_message <<                                                        \
          "', expected '" <<                                                                \
          (message) << "'\n";                                                                 \
        TEST::failed_lines_list.push_back(TEST::test_line);                               \
        break;                                                                            \
      }                                                                                   \
    }                                                                                     \
  }

/** @brief Create a temporary filename.

 This macro assigns a new temporary filename to the string variable given as
 its argument. The filename is created using the filename of the test and the
 line number where this macro is invoked, for example 'Matrix_test.cpp' might
 create a temporary file 'Matrix_test_268.tmp' if NEW_TMP_FILE is used in
 line 268.  All temporary files are deleted if #END_TEST is called.  @p
 filename string will contain the filename on completion of the macro.

 There is a version that defines the extension and one that uses tmp.

 @hideinitializer
 */
#define NEW_TMP_FILE_EXT(filename, extension) filename = TEST::createTmpFileName(__FILE__, __LINE__, extension);


#define NEW_TMP_FILE(filename) filename = TEST::createTmpFileName(__FILE__, __LINE__);


/** @brief Skip the remainder of the current subtest.

 If the condition is not fulfilled, the remainder of the current subtest is
 skipped over. The TEST status is set to FAIL.

 @hideinitializer
 */
#define ABORT_IF(condition)                                                               \
  if (condition)                                                                          \
  {                                                                                       \
    {                                                                                     \
      TEST::test_line = __LINE__;                                                         \
      TEST::this_test = false;                                                            \
      TEST::test = TEST::test && TEST::this_test;                                         \
      TEST::failed_lines_list.push_back(TEST::test_line);                                 \
      TEST::initialNewline();                                                             \
      stdcout << " -  line " << TEST::test_line <<                                        \
        ":  ABORT_IF(" # condition "):  TEST ABORTED\n";                                  \
    }                                                                                     \
    break;                                                                                \
  }

/**
 @brief Print a status message.

 If tests require longer preparations, #STATUS may be used to print some
 intermediate progress messages.  #STATUS uses <code>cout</code> to print
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
#define STATUS(message)                                                                   \
  {                                                                                       \
    TEST::initialNewline();                                                               \
    stdcout << "    line "                                                                \
              <<  __LINE__                                                                \
              << ": status:  "                                                            \
              << message                                                                  \
              << "\n";                                                                    \
  }

/**
 @brief Sets an additional text that is displayed after final result of the test.

 This can be used to provide additional information about the test to the user.
 It is e.g. used to indicate that the DB test were skipped, when there are no
 credentials given.

 @hideinitializer
 */
#define ADD_MESSAGE(message)                                                              \
  TEST::add_message = message;

/**
 @brief Macro that suppresses the warning issued when no subtests are performed

 Please use this macro only if the method cannot be tested at all or cannot be
 tested properly on its own. In the later case, the method must however be tested
 in tests of related methods.  See also @em test_count.

 @hideinitializer
 */
#define NOT_TESTABLE                                                                      \
  TEST::test_count = 1;

//@} // end of ClassTest
