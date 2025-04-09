// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl, Chris Bielow, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <boost/math/special_functions/fpclassify.hpp>

#include <iomanip>
#include <fstream>

#include <QtCore/QFileInfo>

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
      std::string exception_message = "";
      std::string exception_name = "";
      std::string fuzzy_message;
      std::string test_name = "";
      std::vector<std::string> tmp_file_list;
      std::vector<UInt> failed_lines_list;
      StringList whitelist;

      void mainInit(const char* version, const char* class_name, int argc, const char* argv0)
      {
        // if env var "OPENMS_TEST_VERBOSE=True" enable output of successfull line
        char* pverbose = std::getenv("OPENMS_TEST_VERBOSE");
        if (pverbose != nullptr)
        {
          if (std::string(pverbose) == "True") TEST::verbose = 2;
        }

        OpenMS::UniqueIdGenerator::setSeed(2453440375);
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
          String TEST_FILE__template_line;
          String TEST_FILE__line;

          while (TEST::infile.good() && TEST::templatefile.good())
          {
            TEST::templatefile.getline(TEST::line_buffer, 65535);
            TEST_FILE__template_line = TEST::line_buffer;
            TEST::infile.getline(TEST::line_buffer, 65535);
            TEST_FILE__line = TEST::line_buffer;
            TEST_FILE__template_line.trim(); // remove leading and trailing whitespaces (ignore CR/LF line endings on Unix)
            TEST_FILE__line.trim();          // remove leading and trailing whitespaces (ignore CR/LF line endings on Unix)
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
        for (OpenMS::Size i = 0; i < TEST::tmp_file_list.size(); ++i)
          {
            if (!OpenMS::File::remove(TEST::tmp_file_list[i]))
            {
              stdcout << "Warning: unable to remove temporary file '"
                      << TEST::tmp_file_list[i]
                      << "'"
                      << std::endl;
            }
          }
      }

      void
      setWhitelist(const char* const /* file */, const int line,
                   const std::string& whitelist_)
      {
        TEST::whitelist = ListUtils::create<String>(whitelist_);

        if ((TEST::verbose > 1) || (!TEST::this_test && (TEST::verbose > 0)))
        {
          TEST::initialNewline();
          stdcout << " +  line " << line << ":  WHITELIST(\"" << whitelist_
                    << "\"):   whitelist is: " << TEST::whitelist << std::endl;
        }
        return;
      }

      void
      initialNewline()
      {
        if (!newline)
        {
          newline = true;
          std::cout << std::endl;
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

      bool
      validate(const std::vector<std::string>& file_names)
      {
        std::cout << "checking (created temporary files)..." << std::endl;
        bool passed_all = true;
        for (Size i = 0; i < file_names.size(); ++i)
        {
          if (File::exists(file_names[i]))
          {
            FileTypes::Type type = FileHandler::getType(file_names[i]);
            bool passed_single = true;
            bool skipped = false;
            switch (type)
            {
            case FileTypes::MZML:
            {
              if (!MzMLFile().isValid(file_names[i]))
              {
                std::cout << " - Error: mzML file does not validate against XML schema '" << file_names[i] << "'" << std::endl;
                passed_single = false;
              }
              StringList errors, warnings;
              if (!MzMLFile().isSemanticallyValid(file_names[i], errors,
                                                  warnings))
              {
                std::cout << " - Error: mzML file semantically invalid '" << file_names[i] << "'" << std::endl;
                for (Size j = 0; j < errors.size(); ++j)
                {
                  std::cout << "Error - " << errors[j] << std::endl;
                }
                passed_single = false;
              }
            }
            break;

            case FileTypes::MZDATA:
              if (!MzDataFile().isValid(file_names[i], std::cerr))
              {
                std::cout << " - Error: Invalid mzData file '" << file_names[i] << "'" << std::endl;
                passed_single = false;
              }
              break;

            case FileTypes::MZXML:
              if (!MzXMLFile().isValid(file_names[i], std::cerr))
              {
                std::cout << " - Error: Invalid mzXML file '" << file_names[i] << "'" << std::endl;
                passed_single = false;
              }
              break;

            case FileTypes::FEATUREXML:
              if (!FeatureXMLFile().isValid(file_names[i], std::cerr))
              {
                std::cout << " - Error: Invalid featureXML file '" << file_names[i] << "'" << std::endl;
                passed_single = false;
              }
              break;

            case FileTypes::IDXML:
              if (!IdXMLFile().isValid(file_names[i], std::cerr))
              {
                std::cout << " - Error: Invalid idXML file '" << file_names[i] << "'" << std::endl;
                passed_single = false;
              }
              break;

            case FileTypes::CONSENSUSXML:
              if (!ConsensusXMLFile().isValid(file_names[i], std::cerr))
              {
                std::cout << " - Error: Invalid consensusXML file '" << file_names[i] << "'" << std::endl;
                passed_single = false;
              }
              break;

            case FileTypes::INI:
              if (!ParamXMLFile().isValid(file_names[i], std::cerr))
              {
                std::cout << " - Error: Invalid Param file '" << file_names[i] << "'" << std::endl;
                passed_single = false;
              }
              break;

            case FileTypes::TRANSFORMATIONXML:
              if (!TransformationXMLFile().isValid(file_names[i], std::cerr))
              {

                passed_single = false;
              }
              break;

            default:
              skipped = true;
              break;
            }
            //output for single file
            if (skipped)
            {
              std::cout << " +  skipped file '" << file_names[i] << "' (type: " << FileTypes::typeToName(type) << ")" << std::endl;
            }
            else if (passed_single)
            {
              std::cout << " +  valid file '" << file_names[i] << "' (type: " << FileTypes::typeToName(type) << ")" << std::endl;
            }
            else
            {
              passed_all = false;
              std::cout << " -  invalid file '" << file_names[i] << "' (type: " << FileTypes::typeToName(type) << ")" << std::endl;
            }
          }
        }
        //output for all files        
        if (passed_all)
        {
          std::cout << ": passed" << std::endl << std::endl;
        }
        else
        {
          std::cout << ": failed" << std::endl << std::endl;
        }
        return passed_all;
      }

      std::string
      createTmpFileName(const std::string& file, int line, const std::string& extension)
      {
        QFileInfo fi(file.c_str());
        std::string filename = (String(fi.baseName())) + '_' + String(line) + ".tmp" + extension;
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
                           long double number_1, const char* number_1_stringified, bool number_1_is_realtype, Int number_1_written_digits,
                           long double number_2, const char* number_2_stringified, bool /* number_2_is_realtype */, Int number_2_written_digits
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
                    << std::endl;
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
                          << std::setprecision(number_2_written_digits) << number_2 << std::endl;
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
                        << "], message: \"" << TEST::fuzzy_message << "\"" << std::endl;
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
                    << "\"" << std::endl;            
            }
          }
          else
          {
            stdcout << " -  line " << line << ":  TEST_STRING_EQUAL("
                      << string_1_stringified << ',' << string_2_stringified
                      << "): got \"" << string_1 << "\", expected \"" << string_2
                      << "\"" << std::endl;
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
        std::exception_ptr ex = std::current_exception();
        try
        {
          std::rethrow_exception(ex);
        } 
        catch (::OpenMS::Exception::BaseException& e)
        {
          TEST::this_test = false;
          TEST::test = false;
          TEST::all_tests = false;
          {
            TEST::initialNewline();
            out << "Error: Caught unexpected OpenMS exception of type '" << e.getName() << "'";
            if ((e.getLine() > 0) && std::strcmp(e.getFile(), ""))
            {
              out << " thrown in line " << e.getLine() << " of file '" << e.getFile() << "' in function '" << e.getFunction() << "'";
            }
            out << " - Message: " << e.what() << std::endl;
          }
        } /* catch std:: exceptions */
        catch (std::exception& e)
        {
          TEST::this_test = false;
          TEST::test = false;
          TEST::all_tests = false;
          {
            TEST::initialNewline();
            out << "Error: Caught unexpected std::exception\n";
            out << " - Message: " << e.what() << std::endl;
          }
        } /* catch all other exceptions */
        catch (...)
        {
          TEST::this_test = false;
          TEST::test = false;
          TEST::all_tests = false;
          {
            TEST::initialNewline();
            out << "Error: Caught unidentified and unexpected exception - No message." << std::endl;
          }
        }
      }

      int endTestPostProcess(std::ostream& out)
      {
        /* check validity of temporary files if known */
        if (!TEST::validate(TEST::tmp_file_list))
        {
          TEST::all_tests = false;
        }
        if (TEST::verbose == 0)
        {
          out << "Output of successful tests were suppressed. Set the environment variable 'OPENMS_TEST_VERBOSE=True' to enable them." << std::endl;
        } /* check for exit code */
        if (!TEST::all_tests)
        {
          out << "FAILED\n";
          if (TEST::add_message != "")
            out << "Message: " << TEST::add_message << '\n';
          out << "Failed lines: ";
          for (OpenMS::Size i = 0; i < TEST::failed_lines_list.size(); ++i)
          {
            out << TEST::failed_lines_list[i] << " ";
          }
          out << std::endl;
          return 1;
        }
        else
        { /* remove temporary files*/
          TEST::removeTempFiles();
          out << "PASSED";
          if (TEST::add_message != "")
            out << " (" << TEST::add_message << ")";
          out << std::endl;
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
          if (OpenMS::String(TEST::test_name).has('~'))
            out << "Warning: no subtests performed in '" << TEST::test_name << "' (line " << line << ")!\n";
        }
        stdcout << std::endl;
      }
  }
