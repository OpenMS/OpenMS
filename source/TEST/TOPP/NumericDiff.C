// -*- Mode: C++; tab-width: 2; -*-
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

#include <cstdlib> // for strtod()
#include <fstream>
#include <iostream>
#include <ctype.h> // for isspace()
#include <limits> // for NaN
#include <sstream>

/// @Todo Think about another distance than the ratio because of infinite ratios for really small values? (e.g. 0 and 10^-15)? (Clemens)

int argc = 0;
char ** argv = 0;

std::ifstream input_1;
std::ifstream input_2;

std::stringstream line_1;
std::stringstream line_2;

std::ios::pos_type line_1_pos;
std::ios::pos_type line_2_pos;

/// Maximum ratio of numbers allowed
double ratio_max_allowed = 1.0;

/// Maximum ratio of numbers observed so far
double ratio_max = 1.0;

/// Maximum absolute difference of numbers allowed
double absdiff_max_allowed = 12345;

/// Maximum difference of numbers observed so far
double absdiff_max = 0.0;

double number_1 = 0;
char letter_1 = 0;
bool is_number_1 = false;
bool is_space_1 = false;

double number_2 = 0;
char letter_2 = 0;
bool is_number_2 = false;
bool is_space_2 = false;

bool is_absdiff_small = false; 

int line_num_1 = 0;
int line_num_2 = 0;

int line_num_1_max = -1;
int line_num_2_max = -1;

int verbose = 2; // default == 2,  -q == 1,  -Q == 0

std::string line_str_1;
std::string line_str_2;

std::string line_str_1_max;
std::string line_str_2_max;

/// Write out usage information.
void usage()
{
  std::cerr <<
    "Usage: " << argv[0] << " file1 file2 ratio absdiff [ -q | -Q ]\n"
    "\n"
    "where\n"
    "\tfile1:      first input file\n"
    "\tfile2:      second input file\n"
    "\tratio:      acceptable relative error (a number >= 1.0)\n"
    "\tabsdiff:    acceptable absolute error (a number >= 0.0)\n"
		"\t-q:         quiet mode (no output unless differences detected)\n"
		"\t-Q:         very quiet mode (absolutely no output)\n"
		"\n"
		"A non-zero exit status indicates that errors were found."
    "\n";
  exit(1);
  return;
}

void report ( char* message = "<no message>" )
{
	if ( verbose >= 1 )
	{
		std::cout <<
			std::boolalpha <<
			"\n"
			"Error: '" << message << "'\n\n"
			"  file:\t\tinput_1\tinput_2\n"
			" --------------------------------\n"
			"  line_num:\t" << line_num_1 << '\t' << line_num_2 << "\n"
			"  is_number:\t" << is_number_1 << '\t' << is_number_2 << "\n"
			"  is_space:\t" << is_space_1 << '\t' << is_space_2 << "\n"
			"  letters:\t\"" << letter_1 << "\"\t\"" << letter_2 << "\"\n"
			"  numbers:\t" << number_1 << '\t' << number_2 << "\n"
			"\n"
			"  ratio_max: " << ratio_max << "\n"
			"  ratio_max_allowed: " << ratio_max_allowed << "\n\n"
			"  absdiff_max: " << absdiff_max << "\n" 
			"  absdiff_max_allowed: " << absdiff_max_allowed <<
			"\n\n"
			"The offending lines (enclosed in ^$) are:\n"
			"\n" <<
			argv[1] << ':' << line_num_1 << ":\n" <<
			'^' << line_1.str() << "$\n"
			"\n" <<
			argv[2] << ':' << line_num_2 << ":\n" <<
			'^' << line_2.str() << "$\n"
			"\n" <<
			std::endl;
	}
	exit(100);
}


/**@brief Simple diff-like application to compare two input files.
	 Numeric differences are tolerated up to a certain ratio.

   <pre>
	 Usage: NumericDiff file1 file2 ratio

	 where
	   file1:      first input file
	   file2:      second input file
	   ratio:      acceptable relative error (a number >= 1.0)
	   absdiff:    acceptable absolute error (a number >= 0.0)

     -q:         quiet mode (no output unless differences detected)
     -Q:         very quiet mode (absolutely no output)

	 A non-zero exit status indicates that errors were found.
   </pre>
*/
int main ( int main_argc, char ** main_argv)
{
	argc = main_argc;
	argv = main_argv;

	switch ( argc )
	{
	case 4:
	case 5:
		break;
	case 6:
		if ( argv[5] == std::string("-q") )
		{
			verbose = 1;
			break;
		}
		if ( argv[5] == std::string("-Q") )
		{
			verbose = 0;
			break;
		}
	default:
		usage();
		return 1;
	}

	input_1.open(argv[1]);
	if ( !input_1 )
  {
    std::cerr << "Error opening first input file '" << argv[1] <<"'.\n";
    return 2;
  }
	input_1.unsetf(std::ios::skipws);

	input_2.open(argv[2]);
	if ( !input_2 )
  {
    std::cerr << "Error opening second input file '" << argv[2] <<"'.\n";
    return 3;
  }
	input_2.unsetf(std::ios::skipws);

  ratio_max_allowed = std::strtod(argv[3], 0); // returns zero upon failure
  if ( ! ratio_max_allowed )
  {
    std::cerr << "Malformed ratio argument: '" << argv[3] << "'\n";
    return 4;
  }

  absdiff_max_allowed = std::strtod(argv[4], 0); // returns zero upon failure
  if ( absdiff_max_allowed == -12345 )
  {
    std::cerr << "Malformed absdiff argument: '" << argv[4] << "'\n";
    return 5;
  }

	if ( ratio_max_allowed < 1.0 ) ratio_max_allowed = 1/ratio_max_allowed;
	if ( absdiff_max_allowed < 0.0 ) absdiff_max_allowed = -absdiff_max_allowed;

	//------------------------------------------------------------
	// main loop

  while ( input_1 || input_2 )
  {

		{ // read the next line in both files, skip empty lines and lines consisting of whitespace only

			for ( line_str_1.clear(); ++line_num_1, std::getline(input_1,line_str_1); )
			{
				if ( line_str_1.empty() ) continue; // shortcut
				std::string::const_iterator iter = line_str_1.begin(); // loop initialization
				for ( ; iter != line_str_1.end() && isspace(*iter); ++iter ) ; // skip over whitespace
				if ( iter != line_str_1.end() ) break; // line is not empty or whitespace only
			}
			line_1.str(line_str_1);
			line_1.seekp(0);
			line_1.clear();
			// std::cout << line_1.str() << std::endl; // debug

			for ( line_str_2.clear(); ++line_num_2, std::getline(input_2,line_str_2); )
			{
				if ( line_str_2.empty() ) continue; // shortcut
				std::string::const_iterator iter = line_str_2.begin(); // loop initialization
				for ( ; iter != line_str_2.end() && isspace(*iter); ++iter ) ; // skip over whitespace
				if ( iter != line_str_2.end() ) break; // line is not empty or whitespace only
			}
			line_2.str(line_str_2);
			line_2.seekp(0);
			line_2.clear();
			// std::cout << line_2.str() << std::endl; // debug

		}

		while ( line_1 || line_2 )
		{
			is_number_1 = false;
			is_number_2 = false;
			is_space_1 = false;
			is_space_2 = false;
			letter_1 = '\0';
			letter_2 = '\0';
			number_1 = std::numeric_limits<double>::quiet_NaN();
			number_2 = std::numeric_limits<double>::quiet_NaN();

			line_1_pos = line_1.tellg();
			line_1 >> number_1;
			is_number_1 = line_1;

			if ( !is_number_1 )
			{
				line_1.clear();
				line_1.seekg(line_1_pos);
				line_1 >> letter_1;
			}

			is_space_1 = isspace(letter_1);
			if ( is_space_1 )
			{
				line_1 >> std::ws;
			}

			line_2_pos = line_2.tellg();
			line_2 >> number_2;
			is_number_2 = line_2;

			if ( !is_number_2 )
			{
				line_2.clear();
				line_2.seekg(line_2_pos);
				line_2 >> letter_2;
			}

			is_space_2 = isspace(letter_2);
			if ( is_space_2 )
			{
				line_2 >> std::ws;
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
					// not an error, if relative error is small.  We do not drop out of
					// the case distinction here because we want to report the relative
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
								report("number_1 is zero, but number_2 is not");
							}
						}
					}
					else
					{ // number_1 is not zero
						if ( !number_2 )
						{
							if ( !is_absdiff_small )
							{
								report("number_1 is not zero, but number_2 is");
							}
						}
						else
						{ // both numbers are not zero
							double ratio = number_1 / number_2;
							if ( ratio < 0 )
							{
								if ( !is_absdiff_small )
								{
									report("numbers have different signs");
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
											report("ratio of numbers is too large");
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
					report("input_1 is a number, but input_2 is not");
				}
			}
			else
			{ // input_1 is not a number
				if ( is_number_2 )
				{
					report("input_1 is not a number, but input_2 is");
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
							report("input_1 is whitespace, but input_2 is not");
						}
					}
					else
					{ // input_1 is not whitespace
						if ( is_space_2 )
						{
							report("input_1 is not whitespace, but input_2 is");
						}
						else
						{ // both inputs are neither numbers nor whitespace, let us compare them as characters
							if ( letter_1 == letter_2 )
							{ // ok, same characters
								continue;
							}
							else
							{
								report("different letters");
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

			verbose = 10000;
			report
				("This cannot happen.  You should never get here ... "
				 "please report this bug along with the data that produced it."
				);
		}

	}

	if ( verbose >= 2 )
	{
		std::cout <<
			argv[0] << ":  Success.\n\n"
			"  ratio_max: " << ratio_max << "\n"
			"  ratio_max_allowed: " << ratio_max_allowed << "\n\n"
			"  absdiff_max: " << absdiff_max << "\n" 
			"  absdiff_max_allowed: " << absdiff_max_allowed << "\n" <<
			std::endl;

		if ( line_num_1_max == -1 && line_num_2_max == -1 )
		{
			std::cout << "No numeric differences were found.\n" << std::endl;
		}
		else
		{
			std::cout <<
				"Maximum relative error was attained at these lines:\n"
				"(lines are enclosed in ^$)\n" <<
				"\n" <<
				argv[1] << ':' << line_num_1_max << ":\n" <<
				"^"<< line_str_1_max << "$\n"
				"\n"
				"  " << argv[2] << ':' << line_num_2_max << ":\n" <<
				"^"<< line_str_2_max << "$\n" <<
				std::endl;
		}

	}

  return 0;
}


