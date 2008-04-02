// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
			return ratio_max_allowed;
		}

		/// Acceptable relative error (a number >= 1.0)
		void setAcceptableRelative(const double rhs)
		{
			this->ratio_max_allowed = rhs;
			if ( ratio_max_allowed < 1.0 ) ratio_max_allowed = 1/ratio_max_allowed;

		}

		/// Acceptable absolute difference (a number >= 0.0)
		const double & getAcceptableAbsolute() const
		{
			return absdiff_max_allowed;
		}

		/// Acceptable absolute difference (a number >= 0.0)
		void setAcceptableAbsolute(const double rhs)
		{
			this->absdiff_max_allowed = rhs;
			if ( absdiff_max_allowed < 0.0 ) absdiff_max_allowed = -absdiff_max_allowed;
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
			return verbose_level;
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
			this->verbose_level = rhs;
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
		*/
		void setLogDestination(std::ostream & rhs)
		{
			this->log_dest_ = &rhs;
		}

		/**@brief Compare two strings.
		
		This compares all lines of the input.

		returns true in case of success
		*/
		bool compare_strings( std::string const & lhs, std::string const & rhs );
		
		/**@brief Compare two streams of input.
		
		This compares all lines of the input.  Intended to be used for file
		streams.

		returns true in case of success
		*/
		bool compare_streams( std::istream & input_1, std::istream & input_2 );
		
		/**@brief Simple diff-like application to compare two input files.
		Numeric differences are tolerated up to a certain ratio or absolute
		difference.

		where
		@param filename_1 first input file
		@param filename_2 second input file
		@return A non-zero exit status indicates that errors were found.  For the meaning of other numbers, see the code.

		@sa acceptable_ratio
		@sa acceptable_absdiff
		@sa verbose_level
		*/
		bool compare_files( const std::string & filename_1, const std::string & filename_2);

	 protected:

		/**@brief Compare two lines of input.
		
		This implements the core functionality.  Intended to be used for a single
		line of input.

		returns true in case of success
		*/
		bool compare_lines_( std::string const & line_str_1,
												 std::string const & line_str_2
											 );

		/// Report good news.
 		void report_success_() const;

		/// Report bad news.
 		void report_failure_( char const * const message ) const throw(AbortComparison);

    /// Log and results output goes here
		std::ostream * log_dest_;

		/// input_1 name
		std::string input_1_name;
		/// input_2 name
		std::string input_2_name;

		std::stringstream line_1;
		std::stringstream line_2;

		std::ios::pos_type line_1_pos;
		std::ios::pos_type line_2_pos;

		/// Maximum ratio of numbers allowed
		double ratio_max_allowed;

		/// Maximum ratio of numbers observed so far
		double ratio_max;

		/// Maximum absolute difference of numbers allowed
		double absdiff_max_allowed;

		/// Maximum difference of numbers observed so far
		double absdiff_max;

		double number_1;
		char letter_1;
		bool is_number_1;
		bool is_space_1;

		double number_2;
		char letter_2;
		bool is_number_2;
		bool is_space_2;

		bool is_absdiff_small; 

		int line_num_1;
		int line_num_2;

		int line_num_1_max;
		int line_num_2_max;

		int verbose_level;
		
		/**@brief Has comparison been sucessful so far?  Note: this flag is
		changed in report_failure_();
		*/
		bool is_status_success;

		std::string line_str_1_max;
		std::string line_str_2_max;

	}; // class FuzzyStringComparator

}//namespace OpenMS

#endif //OPENMS_CONCEPT_FUZZYSTRINGCOMPARATOR_H

