// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

int argc = 0;
char ** argv = 0;

std::ifstream input_1;
std::ifstream input_2;

std::stringstream line_1;
std::stringstream line_2;

/// Maximum ratio of numbers allowed
double ratio_max_allowed = 0;

/// Maximum ratio of numbers observed so far
double ratio_max = 0;

double number_1 = 0;
char letter_1 = 0;
bool is_number_1 = 0;
bool is_space_1 = 0;

double number_2 = 0;
char letter_2 = 0;
bool is_number_2 = 0;
bool is_space_2 = 0;

int line_num_1 = 0;
int line_num_2 = 0;

int line_num_1_max = 0;
int line_num_2_max = 0;

std::string tmp;

/// Write out usage information.
void usage()
{
  std::cerr <<
    "Usage: " << argv[0] << " file1 file2 ratio \n"
    "\n"
    "where\n"
    "\tfile1:      first input file\n"
    "\tfile2:      second input file\n"
    "\tratio:      acceptable relative error (a number >= 1.0)\n"
    "\n";
  exit(1);
  return;
}

void report ( char* message = "<no message>" )
{
	std::cout <<
		std::boolalpha <<
		"\n"
		"Error: '" << message << "'\n"
		"  line_num:\t" << line_num_1 << '\t' << line_num_2 << "\n"
		"  is_number:\t" << is_number_1 << '\t' << is_number_2 << "\n"
		"  is_space:\t" << is_space_1 << '\t' << is_space_2 << "\n"
		"  letters:\t\"" << letter_1 << "\"\t\"" << letter_2 << "\"\n"
		"  numbers:\t" << number_1 << '\t' << number_2 << "\n"
		"  ratio_max: " << ratio_max << "  ratio_max_allowed: " << ratio_max_allowed <<
		"\n\n" <<
		argv[1] << ':' << line_num_1 << ":\n" <<
		'^' << line_1.str() << "$\n"
		"\n" <<
		argv[2] << ':' << line_num_2 << ":\n" <<
		'^' << line_2.str() << "$\n"
		"\n" <<
		std::endl;
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
   </pre>
*/
int main ( int main_argc, char ** main_argv)
{
	argc = main_argc;
	argv = main_argv;

  if ( argc != 4 ) usage();

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
    std::cerr << "Malformed precision argument: '" << argv[3] << "'\n";
    return 4;
  }

	if ( ratio_max_allowed < 1.0 ) ratio_max_allowed = 1/ratio_max_allowed;

	ratio_max = 0;

	//------------------------------------------------------------
	// main loop

  while ( input_1 || input_2 )
  {

		{ // there should be a better way to do this

			for ( tmp.clear(); ++line_num_1, std::getline(input_1,tmp) && tmp.empty(); ) ;
			line_1.str(tmp);
			line_1.seekp(0);
			line_1.clear();
			// std::cout << line_1.str() << std::endl;

			for ( tmp.clear(); ++line_num_2, std::getline(input_2,tmp) && tmp.empty(); ) ;
			line_2.str(tmp);
			line_2.seekp(0);
			line_2.clear();
			// std::cout << line_2.str() << std::endl;

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

			line_1 >> number_1;
			is_number_1 = line_1;

			if ( !is_number_1 )
			{
				line_1.clear();
				line_1 >> letter_1;
			}

			is_space_1 = isspace(letter_1);
			if ( is_space_1 )
			{
				line_1 >> std::ws;
			}

			line_2 >> number_2;
			is_number_2 = line_2;

			if ( !is_number_2 )
			{
				line_2.clear();
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
					if ( !number_1 )
					{ // number_1 is zero
						if (!number_2 )
						{ // both numbers are zero
							continue;
						}
						else
						{
							report("number_1 is zero, but number_2 is not");
						}
					}
					else
					{ // number_1 is not zero
						if ( !number_2 )
						{
							report("number_1 is not zero, but number_2 is");
						}
						else
						{ // both numbers are not zero
							double ratio = number_1 / number_2;
							if ( ratio < 0 )
							{
								report("numbers have different signs");
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
									if ( ratio_max > ratio_max_allowed )
									{
										report("ratio of numbers is too large");
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
			report
				("This cannot happen.  You should never get here ... "
				 "please report this bug along with the data that produced it."
				);
		}

	}

	std::cout <<
		argv[0] << ":  Success.\n"
		"  ratio_max_allowed: " << ratio_max_allowed << "\n" 
		"  ratio_max: " << ratio_max << "\n"
		"Maximum relative error was attained at these lines:\n" << 
		"  " << argv[1] << ':' << line_num_1_max << "\n"
		"  " << argv[2] << ':' << line_num_2_max << '\n' <<
		std::endl;

  return 0;
}
