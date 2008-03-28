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

#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <sstream>

namespace OpenMS
{
  
  FuzzyStringComparator::FuzzyStringComparator()
    :
    log_dest_(&std::cout),
    input_1_name("input_1"),
    input_2_name("input_2"),
    line_1(),
    line_2(),
    // Maximum ratio of numbers allowed
    ratio_max_allowed(1.0),
    // Maximum ratio of numbers observed so far
    ratio_max(1.0),
    // Maximum absolute difference of numbers allowed
    absdiff_max_allowed(0.0),
    // Maximum difference of numbers observed so far
    absdiff_max(0.0),
    number_1(0),
    letter_1(0),
    is_number_1(false),
    is_space_1(false),
    number_2(0),
    letter_2(0),
    is_number_2(false),
    is_space_2(false),
    is_absdiff_small(false), 
    line_num_1(0),
    line_num_2(0),
    line_num_1_max(-1),
    line_num_2_max(-1),
    // default == 2,  -q == 1,  -Q == 0,  continue after failure == 3
    verbose_level(2),
    is_status_success(true),
    line_str_1_max(),
    line_str_2_max()
  {
  }

  FuzzyStringComparator::~FuzzyStringComparator(){}
  
  void FuzzyStringComparator::report_failure_( char const * const message ) const
    throw (FuzzyStringComparator::AbortComparison)
  {
    // We neither want this entire method be non-const nor make
    // is_status_success a mutable.  So lets hack around it.  (Documented in
    // class.)
    const_cast<bool&>(is_status_success) = false;

    if ( verbose_level >= 1 )
    {
      *log_dest_ <<
				std::boolalpha <<
				"FAILED: '" << message <<
				"'\n\n"
				"  input:\tlhs\trhs\n"
				"  line_num:\t" << line_num_1 << '\t' << line_num_2 << "\n"
				"  col_num:\t" << line_1_pos << '\t' << line_2_pos << "\n"
				" --------------------------------\n"
				"  is_number:\t" << is_number_1 << '\t' << is_number_2 << "\n"
				"  numbers:\t" << number_1 << '\t' << number_2 << "\n"
				"  is_space:\t" << is_space_1 << '\t' << is_space_2 << "\n"
				"  is_letter:\t" << (!is_number_1&&!is_space_1) << '\t' << (!is_number_2&&!is_space_2) << "\n"
				"  letters:\t\"" << letter_1 << "\"\t\"" << letter_2 << "\"\n"
				"  char_codes:\t" << static_cast<UInt>(letter_1) << "\t" << static_cast<UInt>(letter_2) << "\n"
				" --------------------------------\n"
				"  relative_max:        " << ratio_max << "\n"
				"  relative_acceptable: " << ratio_max_allowed << "\n"
				" --------------------------------\n"
				"  absolute_max:        " << absdiff_max << "\n" 
				"  absolute_acceptable: " << absdiff_max_allowed <<
				"\n\n"
				"Offending lines:\n"
				"\n" <<
				input_1_name << ':' << line_num_1 << ':' << line_1_pos << ":\n";
			OpenMS::String pre1(line_1.str());
			pre1 = pre1.prefix(size_t(line_1_pos));
			*log_dest_ << pre1;
			for ( String::iterator iter = pre1.begin(); iter != pre1.end(); ++iter )
				if ( *iter != '\t' ) *iter = ' ';
			*log_dest_ << "!\n" <<
				pre1 << OpenMS::String(line_1.str()).suffix(line_1.str().size()-pre1.size()) << "\n\n" <<
				input_2_name << ':' << line_num_2 << ':' << line_2_pos << ":\n";
			OpenMS::String pre2(line_2.str());
			pre2 = pre2.prefix(size_t(line_2_pos));
			*log_dest_ << pre2;
			for ( String::iterator iter = pre2.begin(); iter != pre2.end(); ++iter )
				if ( *iter != '\t' ) *iter = ' ';
			*log_dest_ << "!\n" <<
				pre2 << OpenMS::String(line_2.str()).suffix(line_2.str().size()-pre2.size()) << "\n\n" <<
				std::endl;
    }

		// If verbose level is low, report only the first error.
    if ( verbose_level < 3 )
    {
      throw FuzzyStringComparator::AbortComparison();
    }
 
    return;
  } // report_failure_()

  void FuzzyStringComparator::report_success_() const
  {
    if (  is_status_success && verbose_level >= 2 )
    {
      *log_dest_ <<
				"PASSED.\n\n"
				"  relative_max:        " << ratio_max << "\n"
				"  relative_acceptable: " << ratio_max_allowed << "\n\n"
				"  absolute_max:        " << absdiff_max << "\n" 
				"  absolute_acceptable: " << absdiff_max_allowed << "\n" <<
				std::endl;

      if ( line_num_1_max == -1 && line_num_2_max == -1 )
      {
				*log_dest_ << "No numeric differences were found.\n" << std::endl;
      }
      else
      {
				*log_dest_ <<
					"Maximum relative error was attained at these lines, enclosed in \"\":\n"
					"\n" <<
					input_1_name << ':' << line_num_1_max << ":\n" <<
					"\""<< line_str_1_max << "\"\n"
					"\n" <<
					input_2_name << ':' << line_num_2_max << ":\n" <<
					"\""<< line_str_2_max << "\"\n" <<
					std::endl;
      }
    } // if verbose_level
    return;
  }
  
  bool FuzzyStringComparator::compare_lines_( std::string const & line_str_1, std::string const & line_str_2 )
  {
		
    line_1.str(line_str_1);
    line_1.seekp(0);
    line_1.clear();
		line_1.unsetf(std::ios::skipws);

    line_2.str(line_str_2);
    line_2.seekp(0);
    line_2.clear();
		line_2.unsetf(std::ios::skipws);

    try
    {
      while ( line_1 && line_2 )
      {
				is_number_1 = false;
				is_number_2 = false;
				is_space_1 = false;
				is_space_2 = false;
				letter_1 = '\0';
				letter_2 = '\0';
				number_1 = std::numeric_limits<double>::quiet_NaN();
				number_2 = std::numeric_limits<double>::quiet_NaN();

				line_1_pos = line_1.tellg(); // save current reading position
				line_1 >> letter_1; // read letter
				// std::cout << ":::" << letter_1 << line_1_pos << std::endl;
				if ( ( is_space_1 = isspace(letter_1) ) ) // is whitespace?
				{
					line_1 >> std::ws; // skip over further whitespace
				}
				else
				{
					line_1.seekg(line_1_pos); // rewind to saved position
					if ( ( is_number_1 = ( line_1 >> number_1 ) ) ) // is a number?
					{
						// letter_1 = '\0';
						// std::cout << line_1_pos << std::endl;
					}
					else
					{
						line_1.clear(); // reset status
						line_1.seekg(line_1_pos); // rewind to saved position
						line_1 >> letter_1; // read letter
					}
				}

				line_2_pos = line_2.tellg(); // save current reading position
				line_2 >> letter_2; // read letter
				if ( ( is_space_2 = isspace(letter_2) ) ) // is whitespace?
				{
					line_2 >> std::ws; // skip over further whitespace
				}
				else
				{
					line_2.seekg(line_2_pos); // rewind to saved position
					if ( ( is_number_2 = ( line_2 >> number_2 ) ) ) // is a number?
					{
						// letter_2 = '\0';
					}
					else
					{
						line_2.clear(); // reset status
						line_2.seekg(line_2_pos); // rewind to saved position
						line_2 >> letter_2; // read letter
					}
				}


				if ( is_number_1 )
				{
					if ( is_number_2 )
					{ // we are comparing numbers

						// check if absolute difference is small
						double absdiff = number_1 - number_2;
						if ( absdiff < 0 )
						{
							absdiff = -absdiff;
						}
						if ( absdiff > absdiff_max )
						{
							absdiff_max = absdiff;
						}
						// If absolute difference is small, large relative errors will be
						// tolerated in the cases below.  But a large absolute difference is
						// not an error, if relative error is small.  We do not jump out of
						// the case distinction here because we want to record the relative
						// error even in case of a successful comparison.
						is_absdiff_small = ( absdiff <= absdiff_max_allowed );
					
						if ( !number_1 )
						{ // number_1 is zero
							if (!number_2 )
							{ // both numbers are zero
								continue;
							}
							else
							{
								if ( !is_absdiff_small )
								{
									report_failure_("number_1 is zero, but number_2 is not");
									continue;
								}
							}
						}
						else
						{ // number_1 is not zero
							if ( !number_2 )
							{
								if ( !is_absdiff_small )
								{
									report_failure_("number_1 is not zero, but number_2 is");
									continue;
								}
							}
							else
							{ // both numbers are not zero
								double ratio = number_1 / number_2;
								if ( ratio < 0 )
								{
									if ( !is_absdiff_small )
									{
										report_failure_("numbers have different signs");
										continue;
									}
								}
								else
								{ // ok, numbers have same sign, but we still need to check their ratio
									if ( ratio < 1 )
									{ // take reciprocal value
										ratio = 1.0 / ratio;
									}
									// by now, we are sure that ratio >= 1
									if ( ratio > ratio_max )
									{ // update running max
										ratio_max = ratio;
										line_num_1_max = line_num_1;
										line_num_2_max = line_num_2;
										line_str_1_max = line_str_1;
										line_str_2_max = line_str_2;
										if ( ratio_max > ratio_max_allowed )
										{
											if ( !is_absdiff_small )
											{
												report_failure_("ratio of numbers is too large");
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
						report_failure_("input_1 is a number, but input_2 is not");
						continue;
					}
				}
				else
				{ // input_1 is not a number
					if ( is_number_2 )
					{
						report_failure_("input_1 is not a number, but input_2 is");
						continue;
					}
					else
					{ // ok, both inputs are not numbers, let us compare them as characters or whitespace
						if ( is_space_1 )
						{
							if ( is_space_2 )
							{ // ok, both inputs are whitespace
								continue;
							}
							else
							{
								report_failure_("input_1 is whitespace, but input_2 is not");
								continue;
							}
						}
						else
						{ // input_1 is not whitespace
							if ( is_space_2 )
							{
								report_failure_("input_1 is not whitespace, but input_2 is");
								continue;
							}
							else
							{ // both inputs are neither numbers nor whitespace, let us compare them as characters
								if ( letter_1 == letter_2 )
								{ // ok, same characters
									continue;
								}
								else
								{
									report_failure_("different letters");
									continue;
								}
							}
						}
					}
				}

				if ( is_absdiff_small )
				{
					is_absdiff_small = false;
					continue;
				}

				verbose_level = 10000;
				report_failure_
					("This cannot happen.  You should never get here ... "
					 "please report this bug along with the data that produced it."
					);

      } // while ( line_1 || line_2 )
			if ( line_1 && !line_2 )
			{
				report_failure_("line from input_2 is shorter than line from input_1");
			}
			if ( !line_1 && line_2 )
			{
				report_failure_("line from input_1 is shorter than line from input_2");
			}
    }
    catch ( FuzzyStringComparator::AbortComparison const& )
    {
      // *log_dest_ << "compare_lines_(): Caught FuzzyStringComparator::AbortComparison\n";
    }

    return is_status_success;
  } // compare_lines_()

  bool FuzzyStringComparator::compare_strings( std::string const & lhs, std::string const & rhs )
  {
		std::istringstream input_1(lhs);
		std::istringstream input_2(rhs);

    std::string line_str_1;
    std::string line_str_2;

    while ( input_1 || input_2 )
    {

			// read the next line in both input streams, skipping over empty lines
			// and lines consisting of whitespace only

			for ( line_str_1.clear(); ++line_num_1, std::getline(input_1,line_str_1); )
			{
				if ( line_str_1.empty() ) continue; // shortcut
				std::string::const_iterator iter = line_str_1.begin(); // loop initialization
				for ( ; iter != line_str_1.end() && isspace(*iter); ++iter ) ; // skip over whitespace
				if ( iter != line_str_1.end() ) break; // line is not empty or whitespace only
			}

			for ( line_str_2.clear(); ++line_num_2, std::getline(input_2,line_str_2); )
			{
				if ( line_str_2.empty() ) continue; // shortcut
				std::string::const_iterator iter = line_str_2.begin(); // loop initialization
				for ( ; iter != line_str_2.end() && isspace(*iter); ++iter ) ; // skip over whitespace
				if ( iter != line_str_2.end() ) break; // line is not empty or whitespace only
			}

			// compare the two lines of input
			if ( !compare_lines_(line_str_1, line_str_2) && verbose_level < 3 ) break;

    } // while ( input_1 || input_2 )
		
		report_success_();
		
		return is_status_success;

  } // compare_strings()

  bool FuzzyStringComparator::compare_streams( std::istream & input_1, std::istream & input_2 )
  {
    std::string line_str_1;
    std::string line_str_2;

    while ( input_1 || input_2 )
    {

			// read the next line in both input streams, skipping over empty lines
			// and lines consisting of whitespace only

			for ( line_str_1.clear(); ++line_num_1, std::getline(input_1,line_str_1); )
			{
				if ( line_str_1.empty() ) continue; // shortcut
				std::string::const_iterator iter = line_str_1.begin(); // loop initialization
				for ( ; iter != line_str_1.end() && isspace(*iter); ++iter ) ; // skip over whitespace
				if ( iter != line_str_1.end() ) break; // line is not empty or whitespace only
			}

			for ( line_str_2.clear(); ++line_num_2, std::getline(input_2,line_str_2); )
			{
				if ( line_str_2.empty() ) continue; // shortcut
				std::string::const_iterator iter = line_str_2.begin(); // loop initialization
				for ( ; iter != line_str_2.end() && isspace(*iter); ++iter ) ; // skip over whitespace
				if ( iter != line_str_2.end() ) break; // line is not empty or whitespace only
			}

			// compare the two lines of input
			if ( !compare_lines_(line_str_1, line_str_2) && verbose_level < 3 ) break;

    } // while ( input_1 || input_2 )

		report_success_();
    
    return is_status_success;

  } // compare_streams()

	bool FuzzyStringComparator::compare_files(const std::string & filename_1, const std::string & filename_2)
  {

		input_1_name = filename_1;
		input_2_name = filename_2;

		if ( input_1_name == input_2_name )
		{
			*log_dest_ << "Error: first and second input file have the same name.  That's cheating!\n";
			return 2;
		}
    
		std::ifstream  input_1_f;
		input_1_f.open(input_1_name.c_str());
		if ( !input_1_f )
		{
			*log_dest_ << "Error opening first input file '" << input_1_name <<"'.\n";
			return 11;
		}
		input_1_f.unsetf(std::ios::skipws);

		std::ifstream  input_2_f;
		input_2_f.open(input_2_name.c_str());
		if ( !input_2_f )
		{
			*log_dest_ << "Error opening second input file '" << input_2_name <<"'.\n";
			return 12;
		}
		input_2_f.unsetf(std::ios::skipws);

		//------------------------------------------------------------
		// main loop

		compare_streams(input_1_f, input_2_f);

		return is_status_success;

	} // compare_files()

} //namespace OpenMS
