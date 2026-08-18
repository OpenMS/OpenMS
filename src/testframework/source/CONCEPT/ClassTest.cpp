// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl, Chris Bielow, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/CONCEPT/FuzzyStringComparator.h>

// The framework depends on the C++ standard library only (see ClassTest.h).
// OpenMS-specific behavior (exception naming, unique-ID seeding) is registered
// by the test projects, see src/tests/class_tests/openms/source/OpenMSTestSupport.cpp.
#include <cstdlib>      // std::getenv, exit
#include <exception>    // std::exception, current-exception handling
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <system_error> // std::error_code

namespace OpenMS::Internal::ClassTest
{
      bool all_tests = true;
      bool equal_files;
      bool newline = false;
      bool test = true;
      bool this_test;
      char line_buffer[65536];
      const char* version_string = nullptr;
      double absdiff = 0.;
      double absdiff_max = 0.;
      double absdiff_max_allowed = 1E-5;
      double ratio = 1.;
      double ratio_max = 1.;
      double ratio_max_allowed = 1. + 1E-5;
      int exception = 0;
      int line_num_1_max = -1;
      int line_num_2_max = -1;
      int start_section_line = 0;
      int test_count = 0;
      int test_line = 0;
      int verbose = 0;
      std::ifstream infile;
      std::ifstream templatefile;
      std::string add_message;
      std::string exception_message;
      std::string exception_name;
      std::string fuzzy_message;
      std::string test_name;
      std::vector<std::string> tmp_file_list;
      std::vector<unsigned int> failed_lines_list;
      std::vector<std::string> whitelist;

      namespace
      {
        // Registered exception translators (see ClassTest.h). A small fixed array is
        // plenty: one per library-under-test project.
        constexpr int max_translators = 8;
        ExceptionTranslator translators[max_translators] = {};
        int translator_count = 0;
      }

      namespace
      {
        bool precondition_tests_enabled = false;
      }

      bool preconditionTestsEnabled()
      {
        return precondition_tests_enabled;
      }

      void setPreconditionTestsEnabled(bool enabled)
      {
        precondition_tests_enabled = enabled;
      }

      void registerExceptionTranslator(ExceptionTranslator translator)
      {
        if (translator == nullptr) return;
        for (int i = 0; i < translator_count; ++i)
        {
          if (translators[i] == translator) return; // idempotent
        }
        if (translator_count < max_translators)
        {
          translators[translator_count++] = translator;
        }
      }

      // NOTE: must only be called while an exception is being handled (i.e. from
      // inside a catch block): both the translators and the fallback below rethrow
      // the current exception with a bare `throw;`, which calls std::terminate if
      // no exception is active.
      std::string describeCaughtException()
      {
        for (int i = 0; i < translator_count; ++i)
        {
          std::string description;
          if (translators[i](description)) return description;
        }
        try
        {
          throw; // rethrow the exception currently being handled
        }
        catch (const std::exception& e)
        {
          return std::string("std::exception - Message: ") + e.what();
        }
        catch (...)
        {
          return "unidentified exception - No message.";
        }
      }

      void mainInit(const char* version, const char* class_name, int argc, const char* argv0)
      {
        // if env var "OPENMS_TEST_VERBOSE=True" enable output of successfull line
        char* pverbose = std::getenv("OPENMS_TEST_VERBOSE");
        if (pverbose != nullptr)
        {
          if (std::string(pverbose) == "True") TEST::verbose = 2;
        }

        TEST::version_string = version;

        if (argc > 1)
        {
          std::cerr
              << "This is " << argv0 << ", the test program for the\n"
              << class_name << " class.\n"
                             "\n"
                             "On successful operation it returns PASSED,\n"
                             "otherwise FAILED is printed.\n";
          exit(1);
        }
      }

      void filesEqual(int line, const char* filename, const char* templatename, const char* filename_stringified, const char* templatename_stringified)
      {
        ++TEST::test_count;
        TEST::test_line = line;

        TEST::equal_files = true;
        TEST::infile.open(filename, std::ios::in);
        TEST::templatefile.open(templatename, std::ios::in);

        if (TEST::infile.good() && TEST::templatefile.good())
        {
          std::string TEST_FILE__template_line;
          std::string TEST_FILE__line;

          while (TEST::infile.good() && TEST::templatefile.good())
          {
            TEST::templatefile.getline(TEST::line_buffer, 65535);
            TEST_FILE__template_line = TEST::line_buffer;
            TEST::infile.getline(TEST::line_buffer, 65535);
            TEST_FILE__line = TEST::line_buffer;
            detail::trim(TEST_FILE__template_line); // remove leading and trailing whitespaces (ignore CR/LF line endings on Unix)
            detail::trim(TEST_FILE__line);          // remove leading and trailing whitespaces (ignore CR/LF line endings on Unix)
            if (TEST_FILE__template_line != TEST_FILE__line)
            {
              TEST::equal_files = false;
              TEST::initialNewline();
              stdcout << "   TEST_FILE_EQUAL: line mismatch:\n    got:      '"
                      << TEST_FILE__line << "'\n    expected: '"
                      << TEST_FILE__template_line << "'\n";
            }
          }
        }
        else
        {
          TEST::equal_files = false;
          {
            TEST::initialNewline();
            stdcout << " +  line "
                    << line
                    << ": TEST_FILE_EQUAL("
                    << filename_stringified
                    << ", "
                    << templatename_stringified;
            stdcout << ") : "
                    << " cannot open file: \"";
            if (!TEST::infile.good())
            {
              stdcout << filename << "\" (input file) ";
            }
            if (!TEST::templatefile.good())
            {
              stdcout << templatename << "\" (template file) ";
            }
            stdcout << "'\n";
          }
        }
        TEST::infile.close();
        TEST::templatefile.close();
        TEST::infile.clear();
        TEST::templatefile.clear();

        TEST::this_test = TEST::equal_files;
        TEST::test = TEST::test && TEST::this_test;
        {
          TEST::initialNewline();
          if (TEST::this_test)
          {
            if (TEST::verbose > 1)
            {
              stdcout << " +  line "
                      << line
                      << ": TEST_FILE_EQUAL("
                      << filename_stringified
                      << ", "
                      << templatename_stringified
                      << "): true";
            }
          }
          else
          {
            stdcout << " -  line "
                    << line
                    << ": TEST_FILE_EQUAL("
                    << filename_stringified
                    << ", "
                    << templatename_stringified
                    << "): false (different files: "
                    << filename
                    << " "
                    << templatename
                    << " )\n";
            TEST::failed_lines_list.push_back(TEST::test_line);
          }
        }
      }

      void removeTempFiles()
      {
        for (std::size_t i = 0; i < TEST::tmp_file_list.size(); ++i)
          {
            // NEW_TMP_FILE only registers the *name*; the test may never have created
            // the file, so a missing file is fine (like OpenMS::File::remove).
            const std::filesystem::path p = detail::to_path(TEST::tmp_file_list[i]);
            std::error_code ec;
            if (std::filesystem::exists(p, ec) && !std::filesystem::remove(p, ec))
            {
              stdcout << "Warning: unable to remove temporary file '"
                      << TEST::tmp_file_list[i]
                      << "'"
                      << '\n';
            }
          }
      }

      void
      setWhitelist(const char* const /* file */, const int line,
                   const std::string& whitelist_)
      {
        TEST::whitelist = detail::split(whitelist_, ',');

        if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))
        {
          TEST::initialNewline();
          stdcout << " +  line " << line << ":  WHITELIST(\"" << whitelist_
                    << "\"):   whitelist is: " << detail::join(TEST::whitelist, ", ") << '\n';
        }
        return;
      }

      void
      initialNewline()
      {
        if (!newline)
        {
          newline = true;
          std::cout << '\n';
        }
        return;
      }

      void
      printWithPrefix(const std::string& text, const int marked)
      {
        std::istringstream is(text);
        std::string line;
        int line_number = 0;
        while (std::getline(is, line))
        {
          ++line_number;
          std::cout << (line_number == marked ? " # :|:  " : "   :|:  ") << line << '\n';
        }
        return;
      }

      std::string
      createTmpFileName(const std::string& file, int line, const std::string& extension)
      {
        std::filesystem::path fp(file);
        std::string filename = (std::string(fp.stem().string())) + '_' + std::to_string(line) + ".tmp" + extension;
        TEST::tmp_file_list.push_back(filename);
        TEST::initialNewline();
        stdcout << "    creating new temporary filename '"
                << filename
                << "' (line "
                << __LINE__
                << ")\n";
        return filename;
      }

      void testRealSimilar(const char* /*file*/, int line,
                           long double number_1, const char* number_1_stringified, bool number_1_is_realtype, int number_1_written_digits,
                           long double number_2, const char* number_2_stringified, bool /* number_2_is_realtype */, int number_2_written_digits
                           )
      {
        TEST::initialNewline();
        ++TEST::test_count;
        TEST::test_line = line;
        TEST::this_test = true;
        if (!number_1_is_realtype)
        {
          TEST::this_test = false;
          stdcout << " -  line " << line << ':'
                    << "TEST_REAL_SIMILAR(" << number_1_stringified << ','
                    << number_2_stringified << "):"
                                     " argument " << number_1_stringified
                    << " does not have a floating point type!  Go fix your code!"
                    << '\n';
          failed_lines_list.push_back(line);
        }
        TEST::test = TEST::test && TEST::this_test;
        if (TEST::this_test)
        {
          TEST::this_test = TEST::isRealSimilar(number_1, number_2);
          TEST::test = TEST::test && TEST::this_test;
          {
            if (TEST::this_test)
            {
              if (TEST::verbose > 1)
              {
                stdcout << " +  line " << line << ":  TEST_REAL_SIMILAR("
                          << number_1_stringified << ',' << number_2_stringified
                          << "): got " << std::setprecision(number_1_written_digits)
                          << number_1 << ", expected "
                          << std::setprecision(number_2_written_digits) << number_2 << '\n';
              }
            }
            else
            {
              stdcout << " -  line " << TEST::test_line
                        << ":  TEST_REAL_SIMILAR(" << number_1_stringified << ','
                        << number_2_stringified << "): got "
                        << std::setprecision(number_1_written_digits) << number_1
                        << ", expected "
                        << std::setprecision(number_2_written_digits) << number_2
                        << " (absolute: " << TEST::absdiff << " ["
                        << TEST::absdiff_max_allowed << "], relative: "
                        << TEST::ratio << " [" << TEST::ratio_max_allowed
                        << "], message: \"" << TEST::fuzzy_message << "\"\n";
              failed_lines_list.push_back(line);
            }
          }
        }
      }

      bool isRealSimilar(long double number_1, long double number_2)
      {
        // Note: The original version of the stuff below was copied from
        // FuzzyStringComparator and then heavily modified for ClassTest.
        // But still the case distinctions should be similar.

        absdiff = 0.;
        ratio = 0.;
        fuzzy_message.clear();

        if (std::isnan(number_1))
        {
          fuzzy_message = "number_1 is nan";
          return false;
        }
        if (std::isnan(number_2))
        {
          fuzzy_message = "number_2 is nan";
          return false;
        }

        // check if absolute difference is small
        absdiff = number_1 - number_2;
        if (absdiff < 0)
        {
          absdiff = -absdiff;
        }
        if (absdiff > absdiff_max)
        {
          absdiff_max = absdiff;
        }
        // If absolute difference is small, large relative errors will be
        // tolerated in the cases below.  But a large absolute difference is
        // not an error, if relative error is small.  We do not jump out of
        // the case distinction here because we want to record the relative
        // error even in case of a successful comparison.
        bool is_absdiff_small = (absdiff <= absdiff_max_allowed);

        if (!number_1) // number_1 is zero
        {
          if (!number_2) // both numbers are zero
          {
            fuzzy_message = "both numbers are zero";
            return true;
          }
          else
          {
            if (!is_absdiff_small)
            {
              fuzzy_message = "number_1 is zero, but number_2 is not small";
              return false;
            }
            else
            {
              fuzzy_message = "number_1 is zero, number_2 is small";
              return true;
            }
          }
        }
        else // number_1 is not zero
        {
          if (!number_2)
          {
            if (!is_absdiff_small)
            {
              fuzzy_message = "number_1 is not zero, but number_2 is";
              return false;
            }
            else
            {
              fuzzy_message = "number_2 is zero, but number_1 is not small";
              return true;
            }
          }
          else // both numbers are not zero
          {
            ratio = number_1 / number_2;
            if (ratio < 0.)
            {
              if (!is_absdiff_small)
              {
                fuzzy_message
                  = "numbers have different signs and difference is not small";
                return false;
              }
              else
              {
                fuzzy_message
                  = "numbers have different signs, but difference is small";
                return true;
              }
            }
            else // ok, numbers have same sign, but we still need to check their ratio
            {
              if (ratio < 1.) // take reciprocal value
              {
                ratio = 1. / ratio;
              }
              // by now, we are sure that ratio >= 1
              if (ratio > ratio_max) // update running max
              {
                ratio_max = ratio;
              }
              if (ratio > ratio_max_allowed)
              {
                if (!is_absdiff_small)
                {
                  fuzzy_message = "ratio of numbers is large";
                  return false;
                }
                else
                {
                  fuzzy_message
                    = "ratio of numbers is large, but numbers are small";
                  return true;
                }
              }
              else
              {
                fuzzy_message = "ratio of numbers is small";
                return true;
              }
            }
          }
        }
      }

      void testStringEqual(const char* /*file*/, int line,
                      const std::string& string_1,
                      const char* string_1_stringified,
                      const std::string& string_2,
                      const char* string_2_stringified)
      {
        ++test_count;
        test_line = line;
        this_test = (string_1 == string_2);
        test = test && this_test;
        {
          initialNewline();
          if (this_test)
          {
            if (TEST::verbose > 1)
            {
            stdcout << " +  line " << line << ":  TEST_STRING_EQUAL("
                    << string_1_stringified << ',' << string_2_stringified
                    << "): got \"" << string_1 << "\", expected \"" << string_2
                    << "\"\n";
            }
          }
          else
          {
            stdcout << " -  line " << line << ":  TEST_STRING_EQUAL("
                      << string_1_stringified << ',' << string_2_stringified
                      << "): got \"" << string_1 << "\", expected \"" << string_2
                      << "\"\n";
            failed_lines_list.push_back(line);
          }
        }
      }

      void testStringSimilar(const char* /*file*/, int line,
                             const std::string& string_1,
                             const char* string_1_stringified,
                             const std::string& string_2,
                             const char* string_2_stringified
                             )
      {
        ++TEST::test_count;
        TEST::test_line = line;

        TEST::fuzzy_message.clear();
        FuzzyStringComparator fsc;
        fsc.setAcceptableAbsolute(absdiff_max_allowed);
        fsc.setAcceptableRelative(ratio_max_allowed);
        fsc.setVerboseLevel(2);
        fsc.setWhitelist(whitelist);
        std::ostringstream os;
        fsc.setLogDestination(os);
        fsc.use_prefix_ = true;

        TEST::this_test = fsc.compareStrings(string_1, string_2);

        TEST::fuzzy_message = os.str();
        TEST::absdiff = fsc.absdiff_max_;
        TEST::ratio = fsc.ratio_max_;
        TEST::line_num_1_max = fsc.line_num_1_max_;
        TEST::line_num_2_max = fsc.line_num_2_max_;

        TEST::test = TEST::test && TEST::this_test;

        TEST::initialNewline();
        if (TEST::this_test)
        {
          if (TEST::verbose > 1)
          {
            stdcout << " +  line " << line << ":  TEST_STRING_SIMILAR("
                      << string_1_stringified << ',' << string_2_stringified << "):  "
                                                                      "absolute: " << TEST::absdiff << " (" << TEST::absdiff_max_allowed
                      << "), relative: " << TEST::ratio << " ("
                      << TEST::ratio_max_allowed << ")    +\n";
            stdcout << "got:\n";
            TEST::printWithPrefix(string_1, TEST::line_num_1_max);
            stdcout << "expected:\n";
            TEST::printWithPrefix(string_2, TEST::line_num_2_max);
          }
        }
        else
        {
          stdcout << " -  line " << TEST::test_line
                    << ": TEST_STRING_SIMILAR(" << string_1_stringified << ','
                    << string_2_stringified << ") ...    -\n"
                                     "got:\n";
          TEST::printWithPrefix(string_1, TEST::line_num_1_max);
          stdcout << "expected:\n";
          TEST::printWithPrefix(string_2, TEST::line_num_2_max);
          stdcout << "message: \n";
          stdcout << TEST::fuzzy_message;
          failed_lines_list.push_back(line);
        }
      }

      bool
      isFileSimilar(const std::string& filename_1,
                    const std::string& filename_2)
      {
        fuzzy_message.clear();
        FuzzyStringComparator fsc;
        fsc.setAcceptableAbsolute(absdiff_max_allowed);
        fsc.setAcceptableRelative(ratio_max_allowed);
        fsc.setVerboseLevel(2);
        fsc.setWhitelist(whitelist);
        std::ostringstream os;
        fsc.setLogDestination(os);
        fsc.use_prefix_ = true;

        bool result = fsc.compareFiles(filename_1, filename_2);

        fuzzy_message = os.str();
        absdiff = fsc.absdiff_max_;
        ratio = fsc.ratio_max_;
        line_num_1_max = fsc.line_num_1_max_;
        line_num_2_max = fsc.line_num_2_max_;

        return result;
      }


      void printLastException(std::ostream& out)
      {
        TEST::this_test = false;
        TEST::test = false;
        TEST::all_tests = false;
        TEST::initialNewline();
        out << "Error: Caught unexpected " << describeCaughtException() << '\n';
      }

      int endTestPostProcess(std::ostream& out)
      {
        if (TEST::verbose == 0)
        {
          out << "Output of successful tests were suppressed. Set the environment variable 'OPENMS_TEST_VERBOSE=True' to enable them.\n";
        } /* check for exit code */
        if (!TEST::all_tests)
        {
          out << "FAILED\n";
          if (TEST::add_message != "")
            out << "Message: " << TEST::add_message << '\n';
          out << "Failed lines: ";
          for (std::size_t i = 0; i < TEST::failed_lines_list.size(); ++i)
          {
            out << TEST::failed_lines_list[i] << " ";
          }
          out << '\n';
          return 1;
        }
        else
        { /* remove temporary files*/
          TEST::removeTempFiles();
          out << "PASSED";
          if (TEST::add_message != "")
            out << " (" << TEST::add_message << ")";
          out << '\n';
          return 0;
        }
      }

      void endSectionPostProcess(std::ostream& out, const int line)
      {
        TEST::all_tests = TEST::all_tests && TEST::test;
        if (TEST::test)
        {
          out << ": passed\n";
        }
        else
        {
          out << ": failed\n";
        }
        if (TEST::test_count == 0)
        {
          if (std::string(TEST::test_name).find('~') != std::string::npos)
            out << "Warning: no subtests performed in '" << TEST::test_name << "' (line " << line << ")!\n";
        }
        stdcout << '\n';
      }
  }
