// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl, Stephan Aiche $
// $Authors: Clemens Groepl, Stephan Aiche $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_FUZZYSTRINGCOMPARATOR_H
#define OPENMS_CONCEPT_FUZZYSTRINGCOMPARATOR_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>

#include <cstdlib> // for strtod()
#include <fstream>
#include <iostream>
#include <cctype> // for isspace()
#include <limits> // for NaN
#include <sstream>
#include <map>

namespace OpenMS
{
  namespace Internal
  {
    namespace ClassTest
    {
      void OPENMS_DLLAPI
      testStringSimilar( const char * file, int line,
                         const std::string & string_1,
                         const char * string_1_stringified,
                         const std::string & string_2,
                         const char * string_2_stringified );
      bool OPENMS_DLLAPI
      isFileSimilar( const std::string &, const std::string & );
    }
  }

  /**
   @brief Fuzzy comparison of strings, tolerates numeric differences.

   */
  class OPENMS_DLLAPI FuzzyStringComparator
  {

      friend void OPENMS_DLLAPI
          Internal::ClassTest::testStringSimilar(
                                                  const char * file,
                                                  int line,
                                                  const std::string & string_1,
                                                  const char * string_1_stringified,
                                                  const std::string & string_2,
                                                  const char * string_2_stringified );
      friend bool OPENMS_DLLAPI
      Internal::ClassTest::isFileSimilar( const std::string &,
                                          const std::string & );

      /// %Internal exception class.
      struct AbortComparison
      {
      };

    public:

      ///@name the fabulous four
      //@{

      /// Constructor
      FuzzyStringComparator();

      /// Destructor
      virtual
      ~FuzzyStringComparator();

      /// Copy constructor intentionally not implemented
      FuzzyStringComparator( const FuzzyStringComparator& rhs );

      /// Assignment operator intentionally not implemented
      FuzzyStringComparator &
      operator=( const FuzzyStringComparator& rhs );

      //@}

      /// Acceptable relative error (a number >= 1.0)
      const double &
      getAcceptableRelative() const
      {
        return ratio_max_allowed_;
      }

      /// Acceptable relative error (a number >= 1.0)
      void
      setAcceptableRelative( const double rhs )
      {
        this->ratio_max_allowed_ = rhs;
        if ( ratio_max_allowed_ < 1.0 ) ratio_max_allowed_ = 1
            / ratio_max_allowed_;

      }

      /// Acceptable absolute difference (a number >= 0.0)
      const double &
      getAcceptableAbsolute() const
      {
        return absdiff_max_allowed_;
      }

      /// Acceptable absolute difference (a number >= 0.0)
      void
      setAcceptableAbsolute( const double rhs )
      {
        this->absdiff_max_allowed_ = rhs;
        if ( absdiff_max_allowed_ < 0.0 ) absdiff_max_allowed_
            = -absdiff_max_allowed_;
      }

      /// White list.  If both lines contain the same element from this list, they are skipped over.
      const StringList &
      getWhitelist() const
      {
        return whitelist_;
      }

      /// White list.  If both lines contain the same element from this list, they are skipped over.
      StringList &
      getWhitelist()
      {
        return whitelist_;
      }

      /// White list.  If both lines contain the same element from this list, they are skipped over.
      void
      setWhitelist( const StringList& rhs )
      {
        whitelist_ = rhs;
      }

      /**@brief verbose level

       - 0 = very quiet mode (absolutely no output)
       - 1 = quiet mode (no output unless differences detected)
       - 2 = default (include summary at end)
       - 3 = continue after errors
       .
       */
      const int &
      getVerboseLevel() const
      {
        return verbose_level_;
      }

      /**@brief verbose level

       - 0 = very quiet mode (absolutely no output)
       - 1 = quiet mode (no output unless differences detected)
       - 2 = default (include summary at end)
       - 3 = continue after errors
       .
       */
      void
      setVerboseLevel( const int rhs )
      {
        this->verbose_level_ = rhs;
      }

      /**@brief get tab width (for column numbers)
       */
      const int &
      getTabWidth() const
      {
        return tab_width_;
      }

      /**@brief set tab width (for column numbers)
       */
      void
      setTabWidth( const int rhs )
      {
        this->tab_width_ = rhs;
      }

      /**@brief get first column (for column numbers)
       */
      const int &
      getFirstColumn() const
      {
        return first_column_;
      }

      /**@brief set first column (for column numbers)
       */
      void
      setFirstColumn( const int rhs )
      {
        this->first_column_ = rhs;
      }

      /**@brief Log output is written to this destination.

       The default is std::cout.  Use std::ostringstream etc. to save the output
       in a string.
       */
      std::ostream &
      getLogDestination() const
      {
        return *log_dest_;
      }

      /**@brief Log output is written to this destination.

       The default is std::cout.  Use std::ostringstream etc. to save the output
       in a string.

       @internal There seems to be an issue with this under Windows, see comment
       in FuzzyStringComparator_test.C

       */
      void
      setLogDestination( std::ostream & rhs )
      {
        this->log_dest_ = &rhs;
      }

      /**@brief Compare two strings.

       This compares all lines of the input.

       returns true in case of no differences found
       */
      bool
      compareStrings( std::string const & lhs, std::string const & rhs );

      /**@brief Compare two streams of input.

       This compares all lines of the input.  Intended to be used for file
       streams.

       returns true in case of no differences found
       */
      bool
      compareStreams( std::istream & input_1, std::istream & input_2 );

      /**@brief Simple diff-like application to compare two input files.
       Numeric differences are tolerated up to a certain ratio or absolute
       difference.

       where
       @param filename_1 first input file
       @param filename_2 second input file
       @return true in case of no differences found

       @sa ratio_max_allowed_
       @sa absdiff_max_allowed_
       @sa verbose_level_
       */
      bool
      compareFiles( const std::string & filename_1,
                    const std::string & filename_2 );

    protected:

      /**@brief Compare two lines of input.

       This implements the core functionality.  Intended to be used for a single
       line of input.

       returns true (non-zero) in case of success
       */
      bool
      compareLines_( std::string const & line_str_1,
                     std::string const & line_str_2 );

      /// Report good news.
      void
      reportSuccess_() const;

      /// Report bad news.
      /// @exception AbortComparison
      void
      reportFailure_( char const * const message ) const;

			/// Write info about hits in the whitelist
			void writeWhitelistCases_(const std::string& prefix) const;

			/// read the next line in input stream, skipping over empty lines
			/// and lines consisting of whitespace only
			void readNextLine_(std::istream& input_stream, std::string& line_string, int& line_number) const;

			/// opens and checks an input file stream std::ifstream
			bool openInputFileStream_(const std::string & filename, std::ifstream& input_stream) const;

      /// Log and results output goes here
      std::ostream * log_dest_;

			/// Name of first input e.g., filename
      std::string input_1_name_;
			/// Name of second input e.g., filename
      std::string input_2_name_;

			/// Stores information about the current input line (i.e., stream for the line and the current position in the stream)
			struct InputLine {
				std::stringstream line_;
				std::ios::pos_type line_position_;

				InputLine()
					: line_()
				{
				}

				/// Initialize the input line to the passed string
				void setToString(const std::string & s)
				{
					line_.str(s);
					line_.seekp(0);
					line_.clear();
					line_.unsetf(std::ios::skipws);

					line_position_ = line_.tellg();
				}

				/// Save current position of the stream
				void updatePosition()
				{
					line_position_ = (Int(line_.tellg()) != -1 ? line_.tellg() : std::ios::pos_type(line_.str().length())); // save current reading position
				}

				/// Resets the stream to the last saved position
				void seekGToSavedPosition()
				{
					line_.clear(); // reset status
					line_.seekg(line_position_); // rewind to saved position
				}

				/**
					Convert to pointer

					The pointer returned is not intended to be referenced, it just indicates success when none of the error flags are set.

					@return	null pointer if either failbit or badbit of the nested std::stringstream is set. A non-null pointer otherwise.
					*/
				operator void * ( ) const
				{
					return line_.operator void *();
				}
			};

			InputLine input_line_1_;
			InputLine input_line_2_;

			int line_num_1_;
			int line_num_2_;

			int line_num_1_max_;
			int line_num_2_max_;

			std::string line_str_1_max_;
			std::string line_str_2_max_;

      /// Maximum ratio of numbers allowed, see @em ratio_max_.
      double ratio_max_allowed_;

      /// Maximum ratio of numbers observed so far, see @em ratio_max_allowed_.
      double ratio_max_;

      /// Maximum absolute difference of numbers allowed, see @em absdiff_max_.
      double absdiff_max_allowed_;

      /// Maximum difference of numbers observed so far, see @em absdiff_max_allowed_.
      double absdiff_max_;

			/// Stores information about characters, numbers, and whitesspaces loaded from the InputStream
			struct StreamElement_ {
				double number;
				unsigned char letter;
				bool is_number;
				bool is_space;

				StreamElement_()
					:	number(0),
						letter(0),
						is_number(false),
						is_space(false)
				{}

				/// reset all elements of the element to default value
				void reset()
				{
					is_number = false;
					is_space = false;
					letter = '\0';
					number = std::numeric_limits<double>::quiet_NaN();
				}


				/// Read the next element from an InputLine and update the InputLine accordingly
				void fillFromInputLine(InputLine& input_line )
				{
					// first reset all internal variables so we do not mess with
					// old values
					reset();

					input_line.updatePosition();
					input_line.line_ >> letter; // read letter
					if ( ( is_space = (isspace(letter)!=0) ) ) // is whitespace?
					{
						input_line.line_ >> std::ws; // skip over further whitespace
					}
					else
					{
						input_line.seekGToSavedPosition();
						if ( ( is_number = ( ( input_line.line_ >> number )!=0) ) ) // is a number?
						{
						}
						else
						{
							input_line.seekGToSavedPosition();
							input_line.line_ >> letter; // read letter
						}
					}
				}
			};

			/// Stores information about characters, numbers, and whitesspaces loaded from the first input stream
			StreamElement_ element_1_;
			/// Stores information about characters, numbers, and whitesspaces loaded from the second input stream
			StreamElement_ element_2_;

			/// Wrapper for the prefix information computed for the failure report
			struct PrefixInfo_ {
				OpenMS::String prefix;
				OpenMS::String prefix_whitespaces;
				int line_column;

				PrefixInfo_(const InputLine& input_line, const int tab_width_, const int first_column_)
					: prefix(input_line.line_.str()), line_column(0)
				{
					prefix = prefix.prefix( size_t(input_line.line_position_) );
					prefix_whitespaces = prefix;
					for ( String::iterator iter = prefix_whitespaces.begin(); iter != prefix_whitespaces.end(); ++iter )
					{
						if ( *iter != '\t' )
						{
							*iter = ' ';
							++line_column;
						}
						else
						{
							line_column = (line_column/tab_width_+1)*tab_width_;
						}
					}
					line_column += first_column_;
				}
			};

      bool is_absdiff_small_;

      int verbose_level_;
      int tab_width_;
      int first_column_;

      /**@brief Has comparison been sucessful so far?  Note: this flag is
       changed in reportFailure_();
       */
      bool is_status_success_;

      /// use a prefix when reporting
      bool use_prefix_;

      StringList whitelist_;
      std::map<String,UInt> whitelist_cases_;

  }; // class FuzzyStringComparator

}//namespace OpenMS

#endif //OPENMS_CONCEPT_FUZZYSTRINGCOMPARATOR_H
