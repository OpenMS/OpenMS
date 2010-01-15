// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
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

      /// Log and results output goes here
      std::ostream * log_dest_;

      /// input_1 name
      std::string input_1_name_;
      /// input_2 name
      std::string input_2_name_;

      std::stringstream line_1_;
      std::stringstream line_2_;

      std::ios::pos_type line_1_pos_;
      std::ios::pos_type line_2_pos_;

      /// Maximum ratio of numbers allowed, see @em ratio_max_.
      double ratio_max_allowed_;

      /// Maximum ratio of numbers observed so far, see @em ratio_max_allowed_.
      double ratio_max_;

      /// Maximum absolute difference of numbers allowed, see @em absdiff_max_.
      double absdiff_max_allowed_;

      /// Maximum difference of numbers observed so far, see @em absdiff_max_allowed_.
      double absdiff_max_;

      double number_1_;
      unsigned char letter_1_;
      bool is_number_1_;
      bool is_space_1_;

      double number_2_;
      unsigned char letter_2_;
      bool is_number_2_;
      bool is_space_2_;

      bool is_absdiff_small_;

      int line_num_1_;
      int line_num_2_;

      int line_num_1_max_;
      int line_num_2_max_;

      int verbose_level_;
      int tab_width_;
      int first_column_;

      /**@brief Has comparison been sucessful so far?  Note: this flag is
       changed in reportFailure_();
       */
      bool is_status_success_;

      std::string line_str_1_max_;
      std::string line_str_2_max_;

      /// use a prefix when reporting
      bool use_prefix_;

      StringList whitelist_;
      std::map<String,UInt> whitelist_cases_;

  }; // class FuzzyStringComparator

}//namespace OpenMS

#endif //OPENMS_CONCEPT_FUZZYSTRINGCOMPARATOR_H
