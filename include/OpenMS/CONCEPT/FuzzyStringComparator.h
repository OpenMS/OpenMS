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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_FUZZYSTRINGCOMPARATOR_H
#define OPENMS_CONCEPT_FUZZYSTRINGCOMPARATOR_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <cstdlib> // for strtod()
#include <fstream>
#include <iostream>
#include <ctype.h> // for isspace()
#include <limits> // for NaN
#include <sstream>

namespace OpenMS
{
	/**
	@brief Fuzzy comparison of strings, tolerates numeric differences.

	*/
	class FuzzyStringComparator
	{
		/// Internal exeption class.
		struct AbortComparison{};
		
	 public:

		///@name the fabulous four
		//@{

		/// Constructor
		FuzzyStringComparator();

		/// Destructor
		virtual ~FuzzyStringComparator();

		/// Copy constructor intentionally not implemented
		FuzzyStringComparator(const FuzzyStringComparator& rhs);

		/// Assignment operator intentionally not implemented
		FuzzyStringComparator & operator=(const FuzzyStringComparator& rhs);

		//@}

		/// Acceptable relative error (a number >= 1.0)
		const double & getAcceptableRelative() const
		{
			return ratio_max__allowed_;
		}

		/// Acceptable relative error (a number >= 1.0)
		void setAcceptableRelative(const double rhs)
		{
			this->ratio_max__allowed_ = rhs;
			if ( ratio_max__allowed_ < 1.0 ) ratio_max__allowed_ = 1/ratio_max__allowed_;

		}

		/// Acceptable absolute difference (a number >= 0.0)
		const double & getAcceptableAbsolute() const
		{
			return absdiff_max__allowed_;
		}

		/// Acceptable absolute difference (a number >= 0.0)
		void setAcceptableAbsolute(const double rhs)
		{
			this->absdiff_max__allowed_ = rhs;
			if ( absdiff_max__allowed_ < 0.0 ) absdiff_max__allowed_ = -absdiff_max__allowed_;
		}

		/**@brief verbose level

		- 0 = very quiet mode (absolutely no output) (was: -Q)
		- 1 = quiet mode (no output unless differences detected) (was: -q)
		- 2 = default (include summary at end)
		- 3 = continue after errors
		.
		*/
		const int & getVerboseLevel() const
		{
			return verbose_level_;
		}

		/**@brief verbose level

		- 0 = very quiet mode (absolutely no output) (was: -Q)
		- 1 = quiet mode (no output unless differences detected) (was: -q)
		- 2 = default (include summary at end)
		- 3 = continue after errors
		.
		*/
		void setVerboseLevel(const int rhs)
		{
			this->verbose_level_ = rhs;
		}

		/**@brief get tab width (for column numbers)
		*/
		const int & getTabWidth() const
		{
			return tab_width_;
		}

		/**@brief set tab width (for column numbers)
		*/
		void setTabWidth(const int rhs)
		{
			this->tab_width_ = rhs;
		}

		/**@brief get first column (for column numbers)
		*/
		const int & getFirstColumn() const
		{
			return first_column_;
		}

		/**@brief set first column (for column numbers)
		*/
		void setFirstColumn(const int rhs)
		{
			this->first_column_ = rhs;
		}

		/**@brief Log output is written to this destination.

		The default is std::cout.  Use std::ostringstream etc. to save the output
		in a string.
		*/
		std::ostream & getLogDestination() const
		{
			return *log_dest_;
		}

		/**@brief Log output is written to this destination.

		The default is std::cout.  Use std::ostringstream etc. to save the output
		in a string.
		
		@todo: default for log destination should be std::cout, but strangely this
		doesn't work with Windows, thus I commented out the first subtest. Maybe a
		problem with static initializers? (see FuzzyStringComparator_test.C:
		TEST_EQUAL(fsc.getLogDestination(),std::cout); ) (Clemens, 2008-03-28)

		*/
		void setLogDestination(std::ostream & rhs)
		{
			this->log_dest_ = &rhs;
		}

		/**@brief Compare two strings.
		
		This compares all lines of the input.

		returns true in case of success
		*/
		bool compareStrings( std::string const & lhs, std::string const & rhs );
		
		/**@brief Compare two streams of input.
		
		This compares all lines of the input.  Intended to be used for file
		streams.

		returns true in case of success
		*/
		bool compareStreams( std::istream & input_1, std::istream & input_2 );
		
		/**@brief Simple diff-like application to compare two input files.
		Numeric differences are tolerated up to a certain ratio or absolute
		difference.

		where
		@param filename_1 first input file
		@param filename_2 second input file
		@return A non-zero exit status indicates that errors were found.  For the meaning of other numbers, see the code.

		@sa acceptable_ratio
		@sa acceptable_absdiff
		@sa verbose_level_
		*/
		bool compare_files( const std::string & filename_1, const std::string & filename_2);

	 protected:

		/**@brief Compare two lines of input.
		
		This implements the core functionality.  Intended to be used for a single
		line of input.

		returns true in case of success
		*/
		bool compareLines_( std::string const & line_str_1,
												 std::string const & line_str_2
											 );

		/// Report good news.
 		void reportSuccess_() const;

		/// Report bad news.
 		void reportFailure_( char const * const message ) const throw(AbortComparison);

    /// Log and results output goes here
		std::ostream * log_dest_;

		/// input_1 name
		std::string input_1_name_;
		/// input_2 name
		std::string input_2_name_;

		std::stringstream line_1_;
		std::stringstream line_2_;

		std::ios::pos_type line_1__pos_;
		std::ios::pos_type line_2__pos_;

		/// Maximum ratio of numbers allowed
		double ratio_max__allowed_;

		/// Maximum ratio of numbers observed so far
		double ratio_max_;

		/// Maximum absolute difference of numbers allowed
		double absdiff_max__allowed_;

		/// Maximum difference of numbers observed so far
		double absdiff_max_;

		double number_1_;
		char letter_1_;
		bool is_number_1_;
		bool is_space_1_;

		double number_2_;
		char letter_2_;
		bool is_number_2_;
		bool is_space_2_;

		bool is_absdiff_small_; 

		int line_num_1_;
		int line_num_2_;

		int line_num_1__max_;
		int line_num_2__max_;

		int verbose_level_;
		int tab_width_;
		int first_column_;
		
		/**@brief Has comparison been sucessful so far?  Note: this flag is
		changed in reportFailure_();
		*/
		bool is_status_success_;

		std::string line_str_1_max_;
		std::string line_str_2_max_;

	}; // class FuzzyStringComparator

}//namespace OpenMS

#endif //OPENMS_CONCEPT_FUZZYSTRINGCOMPARATOR_H

