// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/SYSTEM/File.h>
#include <QDir>
#include <sstream>
#include <cstdlib> // for strtod()
#include <cctype> // for isspace()
#include <limits> // for NaN
#include <fstream>
#include <iostream>
#include <iomanip>

namespace OpenMS
{

  FuzzyStringComparator::FuzzyStringComparator() :
    log_dest_(&std::cout),
    input_1_name_("input_1"),
    input_2_name_("input_2"),
    input_line_1_(),
    input_line_2_(),
    line_num_1_(0),
    line_num_2_(0),
    line_num_1_max_(-1),
    line_num_2_max_(-1),
    line_str_1_max_(),
    line_str_2_max_(),
    ratio_max_allowed_(1.0),
    ratio_max_(1.0),
    absdiff_max_allowed_(0.0),
    absdiff_max_(0.0),
    element_1_(),
    element_2_(),
    is_absdiff_small_(false),
    verbose_level_(2),
    tab_width_(8),
    first_column_(1),
    is_status_success_(true),
    use_prefix_(false),
    whitelist_(),
    whitelist_cases_()
  {
  }

  FuzzyStringComparator::~FuzzyStringComparator()
  {
  }

  const double& FuzzyStringComparator::getAcceptableRelative() const
  {
    return ratio_max_allowed_;
  }

  void FuzzyStringComparator::setAcceptableRelative(const double rhs)
  {
    this->ratio_max_allowed_ = rhs;
    if (ratio_max_allowed_ < 1.0)
      ratio_max_allowed_ = 1
                           / ratio_max_allowed_;

  }

  const double& FuzzyStringComparator::getAcceptableAbsolute() const
  {
    return absdiff_max_allowed_;
  }

  void FuzzyStringComparator::setAcceptableAbsolute(const double rhs)
  {
    this->absdiff_max_allowed_ = rhs;
    if (absdiff_max_allowed_ < 0.0)
      absdiff_max_allowed_
        = -absdiff_max_allowed_;
  }

  const StringList& FuzzyStringComparator::getWhitelist() const
  {
    return whitelist_;
  }

  StringList& FuzzyStringComparator::getWhitelist()
  {
    return whitelist_;
  }

  void FuzzyStringComparator::setWhitelist(const StringList& rhs)
  {
    whitelist_ = rhs;
  }

  void FuzzyStringComparator::setMatchedWhitelist(const std::vector< std::pair<std::string, std::string> >& rhs)
  {
    matched_whitelist_ = rhs;
  }

  const std::vector< std::pair<std::string, std::string> >& FuzzyStringComparator::getMatchedWhitelist() const
  {
    return matched_whitelist_;
  }

  const int& FuzzyStringComparator::getVerboseLevel() const
  {
    return verbose_level_;
  }

  void FuzzyStringComparator::setVerboseLevel(const int rhs)
  {
    this->verbose_level_ = rhs;
  }

  const int& FuzzyStringComparator::getTabWidth() const
  {
    return tab_width_;
  }

  void FuzzyStringComparator::setTabWidth(const int rhs)
  {
    this->tab_width_ = rhs;
  }

  const int& FuzzyStringComparator::getFirstColumn() const
  {
    return first_column_;
  }

  void FuzzyStringComparator::setFirstColumn(const int rhs)
  {
    this->first_column_ = rhs;
  }

  std::ostream& FuzzyStringComparator::getLogDestination() const
  {
    return *log_dest_;
  }

  void FuzzyStringComparator::setLogDestination(std::ostream& rhs)
  {
    this->log_dest_ = &rhs;
  }

  void FuzzyStringComparator::reportFailure_(char const* const message) const
  {
    // We neither want this entire method be non-const nor make
    // is_status_success_ a mutable.  So lets hack around it.  (Documented in
    // class.)
    const_cast<bool&>(is_status_success_) = false;

    if (verbose_level_ >= 1)
    {
      PrefixInfo_ prefix1(input_line_1_, tab_width_, first_column_);
      PrefixInfo_ prefix2(input_line_2_, tab_width_, first_column_);

      std::string prefix;
      if (use_prefix_)
      {
        prefix = "   :|:  ";
      }

      *log_dest_ << std::boolalpha <<
        prefix << "FAILED: '" << message << "'\n" <<
        prefix << "\n" <<
        prefix << "  input:\tin1\tin2\n" <<
        prefix << "  line:\t" << line_num_1_ << '\t' << line_num_2_ << "\n" <<
        prefix << "  pos/col:\t" << input_line_1_.line_position_ << '/' << prefix1.line_column << '\t' << input_line_2_.line_position_ << '/' << prefix2.line_column << "\n" <<
        prefix << " --------------------------------\n" <<
        prefix << "  is_number:\t" << element_1_.is_number << '\t' << element_2_.is_number << "\n" <<
        prefix << "  numbers:\t" << element_1_.number << '\t' << element_2_.number << "\n" <<
        prefix << "  is_space:\t" << element_1_.is_space << '\t' << element_2_.is_space << "\n" <<
        prefix << "  is_letter:\t" << (!element_1_.is_number && !element_1_.is_space) << '\t' << (!element_2_.is_number && !element_2_.is_space) << "\n" <<
        prefix << "  letters:\t\"" << element_1_.letter << "\"\t\"" << element_2_.letter << "\"\n" <<
        prefix << "  char_codes:\t" << static_cast<UInt>(element_1_.letter) << "\t" << static_cast<UInt>(element_2_.letter) << "\n" <<
        prefix << " --------------------------------\n" <<
        prefix << "  relative_max:        " << ratio_max_ << "\n" <<
        prefix << "  relative_acceptable: " << ratio_max_allowed_ << "\n" <<
        prefix << " --------------------------------\n" <<
        prefix << "  absolute_max:        " << absdiff_max_ << "\n" <<
        prefix << "  absolute_acceptable: " << absdiff_max_allowed_ << std::endl;

      writeWhitelistCases_(prefix);

      *log_dest_
        << prefix << "\n"
        << prefix << "Offending lines:\t\t\t(tab_width = " << tab_width_ << ", first_column = " << first_column_ << ")\n"
        << prefix << "\n"
        << prefix << "in1:  " << QDir::toNativeSeparators(File::absolutePath(input_1_name_).toQString()).toStdString() << "   (line: " << line_num_1_ << ", position/column: " << input_line_1_.line_position_ << '/' << prefix1.line_column << ")\n"
        << prefix << prefix1.prefix << "!\n"
        << prefix << prefix1.prefix_whitespaces << OpenMS::String(input_line_1_.line_.str()).suffix(input_line_1_.line_.str().size() - prefix1.prefix.size()) << "\n"
        << prefix <<  "\n"
        << prefix << "in2:  " << QDir::toNativeSeparators(File::absolutePath(input_2_name_).toQString()).toStdString() << "   (line: " << line_num_2_ << ", position/column: " << input_line_2_.line_position_ << '/' << prefix2.line_column << ")\n"
        << prefix << prefix2.prefix << "!\n"
        << prefix << prefix2.prefix_whitespaces << OpenMS::String(input_line_2_.line_.str()).suffix(input_line_2_.line_.str().size() - prefix2.prefix.size()) << "\n"
        << prefix << "\n\n"
        << "Easy Access:" << "\n"
        << QDir::toNativeSeparators(File::absolutePath(input_1_name_).toQString()).toStdString() << ':' << line_num_1_ << ":" << prefix1.line_column << ":\n"
        << QDir::toNativeSeparators(File::absolutePath(input_2_name_).toQString()).toStdString() << ':' << line_num_2_ << ":" << prefix2.line_column << ":\n"
        << "\n"
        #ifdef WIN32
        << "TortoiseGitMerge"
        << " /base:\"" << QDir::toNativeSeparators(File::absolutePath(input_1_name_).toQString()).toStdString() << "\""
        << " /mine:\"" << QDir::toNativeSeparators(File::absolutePath(input_2_name_).toQString()).toStdString() << "\""
        #else
        << "diff"
        << " " << QDir::toNativeSeparators(File::absolutePath(input_1_name_).toQString()).toStdString()
        << " " << QDir::toNativeSeparators(File::absolutePath(input_2_name_).toQString()).toStdString()
        #endif
        << std::endl;
    }

    // If verbose level is low, report only the first error.
    if (verbose_level_ < 3)
    {
      throw FuzzyStringComparator::AbortComparison();
    }

    return;
  } // reportFailure_()

  void FuzzyStringComparator::reportSuccess_() const
  {
    if (is_status_success_ && verbose_level_ >= 2)
    {
      std::string prefix;
      if (use_prefix_)
      {
        prefix = "   :|:  ";
      }

      *log_dest_ <<
        prefix << "PASSED.\n" <<
        prefix << '\n' <<
        prefix << "  relative_max:        " << ratio_max_ << '\n' <<
        prefix << "  relative_acceptable: " << ratio_max_allowed_ << '\n' <<
        prefix << '\n' <<
        prefix << "  absolute_max:        " << absdiff_max_ << '\n' <<
        prefix << "  absolute_acceptable: " << absdiff_max_allowed_ << std::endl;

      writeWhitelistCases_(prefix);

      *log_dest_ << prefix << std::endl;

      if (line_num_1_max_ == -1 && line_num_2_max_ == -1)
      {
        *log_dest_ <<
          prefix << "No numeric differences were found.\n" <<
          prefix << std::endl;
      }
      else
      {
        *log_dest_ <<
          prefix << "Maximum relative error was attained at these lines, enclosed in \"\":\n" <<
          prefix << '\n' <<
          QDir::toNativeSeparators(input_1_name_.c_str()).toStdString() << ':' << line_num_1_max_ << ":\n" <<
          "\"" << line_str_1_max_ << "\"\n" <<
          '\n' <<
          QDir::toNativeSeparators(input_2_name_.c_str()).toStdString() << ':' << line_num_2_max_ << ":\n" <<
          "\"" << line_str_2_max_ << "\"\n" <<
          std::endl;
      }
    }
    return;
  }

  bool FuzzyStringComparator::compareLines_(std::string const& line_str_1, std::string const& line_str_2)
  {

    for (StringList::const_iterator slit = whitelist_.begin();
         slit != whitelist_.end();
         ++slit
         )
    {
      if (line_str_1.find(*slit) != String::npos &&
          line_str_2.find(*slit) != String::npos
          )
      {
        ++whitelist_cases_[*slit];
        // *log_dest_ << "whitelist_ case: " << *slit << '\n';
        return is_status_success_;
      }
    }

    // check matched whitelist
    // If file 1 contains element 1 and file 2 contains element 2, they are skipped over.
    for (std::vector< std::pair<std::string, std::string> >::const_iterator pair_it = matched_whitelist_.begin(); 
         pair_it != matched_whitelist_.end(); ++pair_it)
    {
      if ((line_str_1.find(pair_it->first) != String::npos &&
           line_str_2.find(pair_it->second) != String::npos
          ) ||
          (line_str_1.find(pair_it->second) != String::npos &&
           line_str_2.find(pair_it->first) != String::npos
          )
         )
      {
        // ++whitelist_cases_[*slit];
        // *log_dest_ << "whitelist_ case: " << *slit << '\n';
        return is_status_success_;
      }
    }

    input_line_1_.setToString(line_str_1);
    input_line_2_.setToString(line_str_2);

    try
    {
      while (input_line_1_.ok() && input_line_2_.ok())
      {
        element_1_.fillFromInputLine(input_line_1_);
        element_2_.fillFromInputLine(input_line_2_);

        if (element_1_.is_number)
        {
          if (element_2_.is_number) // we are comparing numbers
          { // check if absolute difference is small
            double absdiff = element_1_.number - element_2_.number;
            if (absdiff < 0)
            {
              absdiff = -absdiff;
            }
            if (absdiff > absdiff_max_)
            {
              absdiff_max_ = absdiff;
            }
            // If absolute difference is small, large relative errors will be
            // tolerated in the cases below.  But a large absolute difference is
            // not an error, if relative error is small.  We do not jump out of
            // the case distinction here because we want to record the relative
            // error even in case of a successful comparison.
            is_absdiff_small_ = (absdiff <= absdiff_max_allowed_);

            if (!element_1_.number) // element_1_.number_ is zero
            {
              if (!element_2_.number) // both numbers are zero
              {
                continue;
              }
              else
              {
                if (!is_absdiff_small_)
                {
                  reportFailure_("element_1_.number_ is zero, but element_2_.number_ is not");
                  continue;
                }
              }
            }
            else // element_1_.number_ is not zero
            {
              if (!element_2_.number)
              {
                if (!is_absdiff_small_)
                {
                  reportFailure_("element_1_.number_ is not zero, but element_2_.number_ is");
                  continue;
                }
              }
              else // both numbers are not zero
              {
                double ratio = element_1_.number / element_2_.number;
                if (ratio < 0)
                {
                  if (!is_absdiff_small_)
                  {
                    reportFailure_("numbers have different signs");
                    continue;
                  }
                }
                else // ok, numbers have same sign, but we still need to check their ratio
                {
                  if (ratio < 1) // take reciprocal value
                  {
                    ratio = 1.0 / ratio;
                  }
                  // by now, we are sure that ratio >= 1
                  if (ratio > ratio_max_) // update running max
                  {
                    ratio_max_ = ratio;
                    line_num_1_max_ = line_num_1_;
                    line_num_2_max_ = line_num_2_;
                    line_str_1_max_ = line_str_1;
                    line_str_2_max_ = line_str_2;
                    if (ratio_max_ > ratio_max_allowed_)
                    {
                      if (!is_absdiff_small_)
                      {
                        reportFailure_("ratio of numbers is too large");
                        continue;
                      }
                    }
                  }
                }
                // okay
                continue;
              }
            }
          }
          else
          {
            reportFailure_("input_1 is a number, but input_2 is not");
            continue;
          }
        }
        else // input_1 is not a number
        {
          if (element_2_.is_number)
          {
            reportFailure_("input_1 is not a number, but input_2 is");
            continue;
          }
          else // ok, both inputs are not numbers, let us compare them as characters or whitespace
          {
            if (element_1_.is_space)
            {
              if (element_2_.is_space) // ok, both inputs are whitespace
              {
                continue;
              }
              else
              {
                if (element_1_.letter == ASCII__CARRIAGE_RETURN) // should be 13 == ascii carriage return char
                {
                  // we skip over '\r'
                  input_line_2_.line_.clear(); // reset status
                  input_line_2_.line_.seekg(input_line_2_.line_position_); // rewind to saved position
                  continue;
                  //reportFailure_("input_1 is carriage return, but input_2_ is not whitespace");
                }
                else
                {
                  reportFailure_("input_1 is whitespace, but input_2 is not");
                }
                continue;
              }
            }
            else // input_1 is not whitespace
            {
              if (element_2_.is_space)
              {
                if (element_2_.letter == ASCII__CARRIAGE_RETURN) // should be 13 == ascii carriage return char
                {
                  // we skip over '\r'
                  input_line_1_.line_.clear(); // reset status
                  input_line_1_.line_.seekg(input_line_1_.line_position_); // rewind to saved position
                  continue;
                  //reportFailure_("input_1 is not whitespace, but input_2 is carriage return");
                }
                else
                {
                  reportFailure_("input_1 is not whitespace, but input_2 is");
                }
                continue;
              }
              else // both inputs are neither numbers nor whitespace, let us compare them as characters
              {
                if (element_1_.letter == element_2_.letter) // ok, same characters
                {
                  continue;
                }
                else
                {
                  reportFailure_("different letters");
                  continue;
                }
              }
            }
          }
        }

        if (is_absdiff_small_)
        {
          is_absdiff_small_ = false;
          continue;
        }

        verbose_level_ = 10000;
        reportFailure_
          ("This cannot happen.  You should never get here ... "
           "please report this bug along with the data that produced it."
          );

      } // while ( input_line_1_ || input_line_2_ )
      if (input_line_1_.ok() && !input_line_2_.ok())
      {
        reportFailure_("line from input_2 is shorter than line from input_1");
      }
      if (!input_line_1_.ok() && input_line_2_.ok())
      {
        reportFailure_("line from input_1 is shorter than line from input_2");
      }
    }
    catch (FuzzyStringComparator::AbortComparison const&)
    {
      // *log_dest_ << "compareLines_(): Caught FuzzyStringComparator::AbortComparison\n";
    }

    return is_status_success_;
  } // compareLines_()

  bool FuzzyStringComparator::compareStrings(std::string const& lhs, std::string const& rhs)
  {
    std::istringstream input_1(lhs);
    std::istringstream input_2(rhs);

    return compareStreams(input_1, input_2);

  } // compareStrings()

  bool FuzzyStringComparator::compareStreams(std::istream& input_1, std::istream& input_2)
  {
    // reset 'success' state to true, in case its currently false due to a prior call (reporting depends on it)
    const_cast<bool&>(is_status_success_) = true;

    std::string line_str_1;
    std::string line_str_2;

    while (input_1 || input_2)
    {

      readNextLine_(input_1, line_str_1, line_num_1_);
      //std::cout << "eof: " << input_1.eof() << " failbit: " << input_1.fail() << " badbit: " << input_1.bad() << " reading " << input_1.tellg () << "chars\n";

      readNextLine_(input_2, line_str_2, line_num_2_);
      //std::cout << "eof: " << input_2.eof() << " failbit: " << input_2.fail() << " badbit: " << input_2.bad() << " reading " << input_2.tellg () << "chars\n";

      // compare the two lines of input
      if (!compareLines_(line_str_1, line_str_2) && verbose_level_ < 3)
        break;

    } // while ( input_1 || input_2 )

    reportSuccess_();

    return is_status_success_;

  } // compareStreams()

  void FuzzyStringComparator::readNextLine_(std::istream& input_stream, std::string& line_string, int& line_number) const
  {
    for (line_string.clear(); static_cast<void>(++line_number), std::getline(input_stream, line_string); )
    {
      if (line_string.empty())
        continue; // shortcut
      std::string::const_iterator iter = line_string.begin(); // loop initialization
      for (; iter != line_string.end() && isspace((unsigned char)*iter); ++iter)
      {
      }
      // skip over whitespace
      if (iter != line_string.end())
        break; // line is not empty or whitespace only
    }
  }

  bool FuzzyStringComparator::compareFiles(const std::string& filename_1, const std::string& filename_2)
  {
    input_1_name_ = filename_1;
    input_2_name_ = filename_2;

    if (input_1_name_ == input_2_name_)
    {
      *log_dest_ << "Error: first and second input file have the same name. That's cheating!\n";
      return false;
    }

    std::ifstream input_1_f;
    if (!openInputFileStream_(input_1_name_, input_1_f))
      return false;

    std::ifstream input_2_f;
    if (!openInputFileStream_(input_2_name_, input_2_f))
      return false;

    //------------------------------------------------------------
    // main loop

    compareStreams(input_1_f, input_2_f);

    return is_status_success_;

  } // compareFiles()

  bool FuzzyStringComparator::openInputFileStream_(const std::string& filename, std::ifstream& input_stream) const
  {
    input_stream.open(filename.c_str(), std::ios::in | std::ios::binary);
    if (!input_stream)
    {
      *log_dest_ << "Error opening first input file '" << filename << "'.\n";
      return false;
    }
    input_stream.unsetf(std::ios::skipws);
    return true;
  }

  void FuzzyStringComparator::writeWhitelistCases_(const std::string& prefix) const
  {
    if (!whitelist_cases_.empty())
    {
      *log_dest_ <<
        prefix << '\n' <<
        prefix << "  whitelist cases:\n";
      Size length = 0;
      for (std::map<String, UInt>::const_iterator wlcit =
             whitelist_cases_.begin(); wlcit != whitelist_cases_.end();
           ++wlcit)
      {
        if (wlcit->first.size() > length)
          length = wlcit->first.size();
      }
      for (std::map<String, UInt>::const_iterator wlcit =
             whitelist_cases_.begin(); wlcit != whitelist_cases_.end();
           ++wlcit)
      {
        *log_dest_ <<
          prefix << "    " << std::setw(length + 3) << std::left <<
          ("\"" + wlcit->first + "\"") << std::setw(3) << std::right <<
          wlcit->second << "x\n";
      }
    }
  }

// we need this only for the stream extraction bug
#ifdef OPENMS_HAS_STREAM_EXTRACTION_BUG
  bool FuzzyStringComparator::StreamElement_::readdigits(InputLine& input_line, std::string& target_buffer, char& c_buffer)
  {
    // we want at least one digit, otherwise it is false
    c_buffer = input_line.line_.peek();
    if (c_buffer == EOF || !isdigit(c_buffer))
    {
      return false;
    }
    else
    {
      target_buffer.push_back(c_buffer);
      c_buffer = input_line.line_.get();
    }

    // try to get more numbers
    while (true)
    {
      c_buffer = input_line.line_.peek();
      if (isdigit(c_buffer) && c_buffer != EOF)
      {
        target_buffer.push_back(c_buffer);
        c_buffer = input_line.line_.get();
      }
      else
      {
        break;
      }
    }

    return true;
  }

  bool FuzzyStringComparator::StreamElement_::tryExtractDouble(InputLine& input_line, double& target)
  {
    char c_peek;
    std::string buffer;

    if (input_line.line_.eof())
      return false;

    c_peek = input_line.line_.peek();

    // .99
    bool have_seen_decimal_separator = false;

    // read sign or first number, otherwise abort
    if ((c_peek == '+' || c_peek == '-' || isdigit(c_peek) || c_peek == '.') && c_peek != EOF)
    {
      buffer.push_back(c_peek);
      c_peek = input_line.line_.get();
      have_seen_decimal_separator = (c_peek == '.');
    }
    else
    {
      return false;
    }

    // try to accumulate more numbers
    readdigits(input_line, buffer, c_peek);

    // if we have no decimal separator, return what we already accumulated
    if (c_peek != '.' && !have_seen_decimal_separator)
    {
      // have we accumulated more then +/-
      if (buffer.size() > 1 || isdigit(buffer[0]))
      {
        target = boost::lexical_cast<double>(buffer);
        return true;
      }
      else
      {
        return false;
      }
    }
    else if (!have_seen_decimal_separator)
    {
      buffer.push_back(c_peek);
      c_peek = input_line.line_.get();

      // try to read the rest of the float
      readdigits(input_line, buffer, c_peek);
    }

    // check if the number is followed by an e+/- something
    if (c_peek == 'e' || c_peek == 'E')
    {
      // add char to buffer
      buffer.push_back(c_peek);
      input_line.line_.get();
      // get next char
      c_peek = input_line.line_.peek();

      // check if we have a sign
      if (c_peek == '+' || c_peek == '-')
      {
        buffer.push_back(c_peek);
        c_peek = input_line.line_.get();
      }

      // try to accumulate the rest of scientific notation
      if (readdigits(input_line, buffer, c_peek))
      {
        target = boost::lexical_cast<double>(buffer);
        return true;
      }
      else
      {
        return false;
      }
    }
    else // just a simple floating point number, convert
    {
      // have we accumulated more then +/-
      if (buffer.size() > 1 || isdigit(buffer[0]))
      {
        target = boost::lexical_cast<double>(buffer);
        return true;
      }
      else
      {
        // not enough content for cast
        return false;
      }
    }
  }

#endif

  void FuzzyStringComparator::StreamElement_::fillFromInputLine(InputLine& input_line)
  {
    // first reset all internal variables so we do not mess with
    // old values
    reset();

    input_line.updatePosition();
    input_line.line_ >> letter; // read letter
    if ((is_space = (isspace(letter) != 0))) // is whitespace?
    {
      input_line.line_ >> std::ws; // skip over further whitespace
    }
    else
    {
      input_line.seekGToSavedPosition();
#ifdef OPENMS_HAS_STREAM_EXTRACTION_BUG
      // is a number? (explicit bool op for C11)
      if (tryExtractDouble(input_line, number))
#else
      if ((is_number = (bool(input_line.line_ >> number))))
#endif
      {
        is_number = true;
      }
      else
      {
        input_line.seekGToSavedPosition();
        input_line.line_ >> letter; // read letter
      }
    }
  }

  void FuzzyStringComparator::StreamElement_::reset()
  {
    is_number = false;
    is_space = false;
    letter = '\0';
    number = std::numeric_limits<double>::quiet_NaN();
  }

  FuzzyStringComparator::StreamElement_::StreamElement_() :
    number(0),
    letter(0),
    is_number(false),
    is_space(false)
  {
  }

  FuzzyStringComparator::InputLine::InputLine() :
    line_()
  {
  }

  void FuzzyStringComparator::InputLine::setToString(const std::string& s)
  {
    line_.str(s);
    line_.seekp(0);
    line_.clear();
    line_.unsetf(std::ios::skipws);

    line_position_ = line_.tellg();
  }

  void FuzzyStringComparator::InputLine::updatePosition()
  {
    line_position_ = (Int(line_.tellg()) != -1 ? line_.tellg() : std::ios::pos_type(line_.str().length())); // save current reading position
  }

  void FuzzyStringComparator::InputLine::seekGToSavedPosition()
  {
    line_.clear(); // reset status
    line_.seekg(line_position_); // rewind to saved position
  }

  bool FuzzyStringComparator::InputLine::ok() const
  {
    return !line_.fail(); // failbit AND badbit are both NOT set; using fail() seems the only portable solution for both C++98 and C++11
                          // operator bool() (C++11 only) and operator void*() (C++98 only) are both not very sexy since they are not "safe bool idiomatic" and would require
                          // a macro here... So we use a real function name (both internally and externally)
  }

  FuzzyStringComparator::PrefixInfo_::PrefixInfo_(const InputLine& input_line, const int this_tab_width_, const int this_first_column_) :
    prefix(input_line.line_.str()), line_column(0)
  {
    prefix = prefix.prefix(size_t(input_line.line_position_));
    prefix_whitespaces = prefix;
    for (String::iterator iter = prefix_whitespaces.begin(); iter != prefix_whitespaces.end(); ++iter)
    {
      if (*iter != '\t')
      {
        *iter = ' ';
        ++line_column;
      }
      else
      {
        line_column = (line_column / this_tab_width_ + 1) * this_tab_width_;
      }
    }
    line_column += this_first_column_;
  }

} //namespace OpenMS
