// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl, Chris Bielow, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

// Std-only: no libOpenMS or third-party headers, so tests of every OpenMS library
// (incl. ones that don't link libOpenMS, like OpenSwathAlgo) can use it and
// recompiles stay cheap. Utilities are its own (ClassTestUtils.h); OpenMS behavior
// is registered by the test projects (openms/source/OpenMSTestSupport.cpp).
#include <OpenMS/CONCEPT/ClassTestUtils.h>
#include <OpenMS/CONCEPT/MacrosTest.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <type_traits>

using XMLCh = char16_t; // Xerces-C++ uses char16_t for UTF-16 strings that we need to output in tests

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

      /// Creates a temporary file name from the test name and the line with the specified extension
      std::string
      createTmpFileName(const std::string& file, int line, const std::string& extension = "");

      // isRealType() and writtenDigits() -- used by #TEST_REAL_SIMILAR -- are provided
      // by ClassTestUtils.h. They classify by convertibility to double, so the framework
      // does not need to know library value types (DataValue, ParamValue, ...) by name.

      /** @brief Compare floating point numbers using @em absdiff_max_allowed and
       @em ratio_max_allowed.

       Side effects: Updates #fuzzy_message.
       */
      void
      testRealSimilar(const char* file, int line, long double number_1,
                      const char* number_1_stringified,
                      bool number_1_is_realtype, int number_1_written_digits,
                      long double number_2, const char* number_2_stringified,
                      bool /* number_2_is_realtype */, int number_2_written_digits);

      /// used by testRealSimilar()
      bool
      isRealSimilar(long double number_1, long double number_2);

      /**@brief Compare strings using @em absdiff_max_allowed and @em ratio_max_allowed.

       This is called by the #TEST_STRING_SIMILAR macro.

       Side effects: Updates #absdiff, #ratio, #fuzzy_message, #line_num_1_max
       and #line_num_2_max.
       */
      void
      testStringSimilar(const char* file, int line,
                        const std::string& string_1,
                        const char* string_1_stringified,
                        const std::string& string_2,
                        const char* string_2_stringified);

      /// used by TEST_STRING_EQUAL
      void
      testStringEqual(const char* file, int line,
                      const std::string& string_1,
                      const char* string_1_stringified,
                      const std::string& string_2,
                      const char* string_2_stringified);

      /**@brief Compare files using @em absdiff_max_allowed and @em ratio_max_allowed.

       Side effects: Updates #absdiff, #ratio, #fuzzy_message, #line_num_1_max
       and #line_num_2_max.
       */
      bool
      isFileSimilar(const std::string& filename_1,
                    const std::string& filename_2);

      /// make sure we have a newline before results from first subtest
      void
      initialNewline();

      /// print the text, each line gets a prefix, the marked line number gets a special prefix
      void
      printWithPrefix(const std::string& text, const int marked = -1);

      /**
         @brief Set up some classtest variables as obtained from the 'START_TEST' macro
                and check that no additional arguments were passed to the test executable.

         @param[in] version A version string, obtained from 'START_TEST(FuzzyStringComparator, "<VERSION>")'
         @param[in] class_name The class under test (used for error messages etc), obtained from 'START_TEST(FuzzyStringComparator, "<VERSION>")'
         @param[in] argc The number of arguments to the main() function of the class test (must be 1; test will quit otherwise)
         @param[in] argv0 Name of the executable (for debug output)
      */
      void mainInit(const char* version, const char* class_name, int argc, const char* argv0);

      /**
        @brief Test if two files are exactly equal (used in TEST_FILE_EQUAL macro)

        @param[in] line The line where the macro was called (for reporting)
        @param[in] filename The temp file
        @param[in] templatename The ground truth file
        @param[in] filename_stringified The expression used as the first macro argument
        @param[in] templatename_stringified The expression used as the second macro argument
      */
      void filesEqual(int line, const char* filename, const char* templatename, const char* filename_stringified, const char* templatename_stringified);

      /// removed all temporary files created with the NEW_TMP_FILE macro
      void removeTempFiles();
      
      /// set the whitelist_
      void
      setWhitelist(const char* const /* file */, const int line,
                   const std::string& whitelist);

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

      /// List of all failed lines for summary at the end of the test
      extern std::vector<unsigned int> failed_lines_list;

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
       #isRealSimilar(), #testStringSimilar(), #isFileSimilar().  Read by
       #TEST_REAL_SIMILAR, #TEST_STRING_SIMILAR, #TEST_FILE_SIMILAR;
       */
      extern std::string fuzzy_message;

      /// (Flags whether a new line is in place, depending on context and verbosity setting.  Used by initialNewline() and some macros.)
      extern bool newline;

      template <typename T1, typename T2>
      void
      testEqual(const char* /*file*/, int line, const T1& expression_1,
                const char* expression_1_stringified,
                const T2& expression_2,
                const char* expression_2_stringified)
      {
        ++test_count;
        test_line = line;
        // For std::string targets, stringify the comparison value rather than
        // constructing std::string(expression_2): the old OpenMS::String had
        // numeric/char converting constructors (e.g. String(114) == "114") that
        // std::string lacks. detail::toString reproduces that behavior (and must
        // format numbers exactly like StringUtils::toStr, see ClassTestUtils.h).
        if constexpr (std::is_same_v<T1, std::string>)
        {
          this_test = bool(expression_1 == detail::toString(expression_2));
        }
        else
        {
          this_test = bool(expression_1 == T1(expression_2));
        }
        test &= this_test;
        {
          initialNewline();
          if (!this_test || verbose > 1)
          {
            stdcout << ' ' << (this_test ? '+' : '-') << "  line " << line << " : TEST_EQUAL(" << expression_1_stringified << ','
                    << expression_2_stringified << "): got '";

            // we can't print wide chars directly using operator<< so we need to test for it
            if constexpr (std::is_same_v<std::remove_cv_t<T1>, XMLCh*> || std::is_same_v<std::remove_cv_t<T2>, XMLCh*>)
            {
              stdcout << (expression_1 == nullptr ? "(null)" : "(XMLCh*)") << "', expected '"
                      << (expression_2 == nullptr ? "(null)" : "(XMLCh*)") << "'\n";
            }
            else if constexpr (std::is_enum_v<T1> && std::is_enum_v<T2>)
            {
              stdcout << static_cast<int>(expression_1) << "', expected '" << static_cast<int>(expression_2) << "'\n";
            }
            else
            {
              detail::printValue(stdcout, expression_1);
              stdcout << "', expected '";
              detail::printValue(stdcout, expression_2);
              stdcout << "'\n";
            }
          }
          if (!this_test)
          {
            failed_lines_list.push_back(line);
          }
        }
      }

      inline void testTrue(const char* /*file*/, int line, const bool expression_1, const char* expression_1_stringified)
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

      inline void testFalse(const char* /*file*/, int line, const bool expression_1, const char* expression_1_stringified)
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
            else
            {
              detail::printValue(stdcout, expression_1);
              stdcout << "', forbidden is '";
              detail::printValue(stdcout, expression_2);
              stdcout << "'\n";
            }
          }
          if (!this_test)
          {
            failed_lines_list.push_back(line);
          }
        }
      }


      /**
        @brief Translator producing a description of the currently handled exception.

        Called with an exception being handled (rethrow it with @c throw; internally).
        If the translator recognizes the type, it writes a one-line description
        (e.g. name, origin, message) into @p description and returns true;
        otherwise it returns false and must leave @p description alone.

        The framework itself knows only @c std::exception. Test projects register
        translators for their library's exception hierarchy -- see
        src/tests/class_tests/openms/source/OpenMSTestSupport.cpp -- which is how
        failure reports show OpenMS exception names without the framework
        depending on libOpenMS. (Same design as Catch2's exception translators.)
      */
      using ExceptionTranslator = bool (*)(std::string& description);

      /// Register an exception translator (idempotent; bounded, first-come first-served).
      void registerExceptionTranslator(ExceptionTranslator translator);

      /// Whether TEST_PRECONDITION/POSTCONDITION_VIOLATED actually run their command.
      /// Off by default; the test-support TU of the library under test enables it when
      /// the library was built with its precondition checks active (OPENMS_ASSERTIONS),
      /// so the expectation always matches the library's real behavior.
      bool preconditionTestsEnabled();

      /// See preconditionTestsEnabled().
      void setPreconditionTestsEnabled(bool enabled);

      /// Describe the exception currently being handled: registered translators first,
      /// then std::exception's what(), then a generic fallback.
      std::string describeCaughtException();

      void printLastException(std::ostream& out);

      int endTestPostProcess(std::ostream& out);

      void endSectionPostProcess(std::ostream& out, const int line);
    }
  }
}

// A namespace alias - apparently these cannot be documented using doxygen (?)
namespace TEST = OpenMS::Internal::ClassTest;

/**
 @defgroup ClassTest Class test macros

 @brief The macros of OpenMS' class-test framework (library
 <code>OpenMSTestFramework</code>), used by the test programs under
 <code>src/tests/class_tests/</code>.

 The framework lives in its own static library, <code>OpenMSTestFramework</code>
 (<code>src/testframework/</code>), and depends on the C++ standard library
 ONLY — it does not link or include libOpenMS. That is what lets the tests of
 every OpenMS library (including ones that do not link libOpenMS, such as
 OpenSwathAlgo) use the same macros, and it keeps the framework working
 unchanged when libOpenMS is split into smaller libraries.

 A test program is one translation unit: #START_TEST / #END_TEST wrap the whole
 program, #START_SECTION / #END_SECTION wrap each subtest, and the elementary
 macros (#TEST_EQUAL, #TEST_REAL_SIMILAR, #TEST_STRING_SIMILAR,
 #TEST_FILE_EQUAL, #TEST_FILE_SIMILAR, #TEST_EXCEPTION, ...) record results.
 On success the program prints "PASSED", otherwise "FAILED" plus the failing
 lines. With <code>-v</code> it prints information about subsections, with
 <code>-V</code> (or the environment variable
 <code>OPENMS_TEST_VERBOSE=True</code>) about every elementary test.

 Because the framework knows nothing about the library under test,
 library-specific behavior is <em>registered</em> by the test project
 (see <code>src/tests/class_tests/openms/source/OpenMSTestSupport.cpp</code>,
 compiled into every openms/openms_gui test executable):

 - OpenMS::Internal::ClassTest::registerExceptionTranslator() — so failure
   reports show the name and origin of unexpected OpenMS exceptions (the
   framework itself only knows <code>std::exception</code>);
 - OpenMS::Internal::ClassTest::setPreconditionTestsEnabled() — so
   #TEST_PRECONDITION_VIOLATED / #TEST_POSTCONDITION_VIOLATED expect a throw
   exactly when the library's #OPENMS_PRECONDITION checks are compiled in;
 - OpenMS::UniqueIdGenerator seeding, so unique IDs written into test output
   files are reproducible against the reference files.

 Values in failure reports print via <code>operator<<</code> found at the
 macro-expansion site (so any type your test can see prints as usual);
 #TEST_REAL_SIMILAR accepts floating point values and class types implicitly
 convertible to <code>double</code> (e.g. DataValue, ParamValue) — no
 registration needed.

 Schema validation of temporary files is explicit: tests that write XML call
 #VALIDATE_TMP_FILES (or #VALIDATE_FILE) from
 <code>OpenMS/TestFileValidation.h</code> in the openms class-test project —
 it needs the FORMAT layer and therefore deliberately does not live in the
 framework.

 To create a test, follow the guidelines in @ref developer_faq (section "How to add a new class test").
 Look at existing test files in src/tests/class_tests/ for examples.

 @ingroup Concept

 */
/** @{ */

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
 any of these exceptions occurs, all tests fail. The report describes the
 caught exception via the registered exception translators (see
 OpenMS::Internal::ClassTest::registerExceptionTranslator), so e.g. %OpenMS
 exceptions are reported with their name and origin.  The #END_TEST macro
 also closes the <code>try</code> block.  This
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
 to mark a subtest that has no matching method declaration.

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
 compile-time errors will occur.  #END_SECTION catches any exception and
 reports it via the registered exception translators (see
 OpenMS::Internal::ClassTest::registerExceptionTranslator).
 After the exception is caught, the execution will continue with
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


/**	@brief Generic equality macro.

 This macro uses the operator == to check its two arguments for equality.
 Besides handling some internal stuff, it basically evaluates ((a) == (b)).

 Remember that operator == has to be defined somehow for the two argument
 types. Additionally the << operator needs to be defined.
 If only == is available you will get a compilation error. As workaround
 use TEST_EQUAL(a==b, true) thereby making bug tracing harder,
 as you won't see the values of a and b.


 @note This macro evaluates its arguments once or twice, depending on verbosity settings.

 @param[in] a value/object to test
 @param[in] b expected value

 @hideinitializer
 */
#define TEST_EQUAL(a, b) TEST::testEqual(__FILE__, __LINE__, (a), (# a), (b), (# b));

/**	@brief Boolean test macro.

 This macro tests if its argument evaluates to 'true'.
 If possible use TEST_EQUAL(a, b) instead of TEST_TRUE(a==b), because the latter makes bug tracing harder.

 @param[in] a value/object convertible to bool
 
 @hideinitializer
*/
#define TEST_TRUE(a) TEST::testTrue(__FILE__, __LINE__, (a), (#a));

/**	@brief Boolean test macro.

 This macro tests if its argument evaluates to 'false'.
 If possible use TEST_NOT_EQUAL(a, b) instead of TEST_FALSE(a!=b), because the latter makes bug tracing harder.

 @param[in] a value/object convertible to bool

 @hideinitializer
*/
#define TEST_FALSE(a) TEST::testFalse(__FILE__, __LINE__, (a), (#a));


/**	@brief Generic inequality macro.

 This macro checks for inequality just like #TEST_EQUAL tests for equality.
 The only difference between the two macros is that #TEST_NOT_EQUAL evaluates
 !((a) == (b)).

 @param[in] a value/object to test
 @param[in] b forbidden value

 @hideinitializer
 */
#define TEST_NOT_EQUAL(a, b) TEST::testNotEqual(__FILE__, __LINE__, (a), (# a), (b), (# b));

/**	@brief std::string equality macro.

 Both arguments are converted to std::string and tested for equality.  (That
 is, we check whether <code>(std::string(a) == std::string(b))</code> holds.)

 @note This macro evaluates its arguments once or twice, depending on verbosity settings.

 @param[in] a value to test
 @param[in] b expected value

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

 @param[in] a value to test
 @param[in] b expected value

 @hideinitializer
 */
#define TEST_REAL_SIMILAR(a, b) TEST::testRealSimilar(__FILE__, __LINE__, (a), (# a), TEST::isRealType(a), TEST::writtenDigits(a), (b), (# b), TEST::isRealType(b), TEST::writtenDigits(b));

/**	@brief std::string similarity macro.

 Compares the two strings using @em FuzzyStringComparator with the settings of
 #TOLERANCE_ABSOLUTE and #TOLERANCE_RELATIVE.

 @note This macro evaluates its arguments once or twice, depending on verbosity settings.

 @note Both arguments are converted to @c std::string.  The actual comparison
 is done by testStringSimilar().

 @param[in] a value to test
 @param[in] b expected value

 @hideinitializer
 */
#define TEST_STRING_SIMILAR(a, b) TEST::testStringSimilar(__FILE__, __LINE__, (a), (# a), (b), (# b));

/**	@brief File similarity macro.

 Compares the two files using @em FuzzyStringComparator with the settings of
 #TOLERANCE_ABSOLUTE and #TOLERANCE_RELATIVE.

 @note This macro evaluates its arguments once or twice, depending on verbosity settings.

 @note The actual comparison is done by isFileSimilar().

 @param[in] a value to test
 @param[in] b expected value

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
                    << TEST::precisionWrapper(TEST::absdiff)                                    \
                    << " ("                                                               \
                    << TEST::precisionWrapper(TEST::absdiff_max_allowed)                        \
                    << "), relative: "                                                    \
                    << TEST::precisionWrapper(TEST::ratio)                                      \
                    << " ("                                                               \
                    << TEST::precisionWrapper(TEST::ratio_max_allowed)                          \
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

 @param[in] exception_type the exception-class
 @param[in] command any general C++ or OpenMS-specific command

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
    catch (...)                                                                           \
    {                                                                                     \
      /* the framework knows only std::exception; registered translators (see  */         \
      /* OpenMSTestSupport.cpp) describe library exception types by name       */         \
      TEST::exception = 2;                                                                \
      TEST::exception_name = TEST::describeCaughtException();                             \
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
      }                                                                                   \
    }                                                                                     \
  }
  
/** @brief Precondition test macro

  This macro checks if a precondition violation is detected while executing the command,
  similar to <code>TEST_EXCEPTION(Exception::Precondition,command)</code>.
  However the test is executed only when the library's #OPENMS_PRECONDITION macros
  are active; the test project registers that fact via setPreconditionTestsEnabled()
  (see OpenMSTestSupport.cpp), keyed on the same OPENMS_ASSERTIONS the library
  uses (see Macros.h).

 @param[in] command any general C++ or OpenMS-specific command

  @hideinitializer
 */
#define TEST_PRECONDITION_VIOLATED(command)                                               \
  {                                                                                       \
    if (TEST::preconditionTestsEnabled())                                                 \
    {                                                                                     \
      TEST_EXCEPTION(Exception::Precondition, command);                                   \
    }                                                                                     \
    else                                                                                  \
    {                                                                                     \
      STATUS("TEST_PRECONDITION_VIOLATED(" # command ")  -  skipped");                    \
    }                                                                                     \
  }

/** @brief Postcondition test macro

  This macro checks if a postcondition violation is detected while executing the command,
  similar to <code>TEST_EXCEPTION(Exception::Postcondition,command)</code>.
  However the test is executed only when the library's #OPENMS_POSTCONDITION macros
  are active; see #TEST_PRECONDITION_VIOLATED for how that is determined.

 @param[in] command any general C++ or OpenMS-specific command

  @hideinitializer
 */
#define TEST_POSTCONDITION_VIOLATED(command)                                               \
  {                                                                                       \
    if (TEST::preconditionTestsEnabled())                                                 \
    {                                                                                     \
      TEST_EXCEPTION(Exception::Postcondition, command);                                   \
    }                                                                                     \
    else                                                                                  \
    {                                                                                     \
      STATUS("TEST_POSTCONDITION_VIOLATED(" # command ")  -  skipped");                    \
    }                                                                                     \
  }


/**	@brief Exception test macro (with test for exception message).

 This macro checks if a given type of exception occurred while executing the
 given command and additionally tests for the message of the exception.

 Example:  #TEST_EXCEPTION_WITH_MESSAGE(Exception::IndexOverflow, vector[-1], "a null pointer was specified")

 If no, a wrong exception occurred or a wrong message is returned, false is
 returned, otherwise true.

 @param[in] exception_type the exception-class
 @param[in] command any general C++ or OpenMS-specific command
 @param[in] message the message the exception should give

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
    catch (...)                                                                           \
    {                                                                                     \
      TEST::exception = 2;                                                                \
      TEST::exception_name = TEST::describeCaughtException();                             \
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

/** @} */ // end of defgroup ClassTest
