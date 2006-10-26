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

#include <cstdlib> // for strtol()
#include <fstream>
#include <iostream>

/// Write out usage information.
void usage( int argc, char ** argv)
{
  std::cerr <<
    "Usage: " << argv[0] << " infile outfile precision\n"
    "\n"
    "where\n"
    "\tinfile:      input file\n"
    "\toutfile:     output file\n"
    "\tprecision:   number of significant digits\n"
    "\n";
  exit(1);
  return;
}

/**@brief Simple minded cat-like application that reads the input file and
   writes it to the output file, but with all numbers rounded to precision
   digits.

   <pre>
   Usage: RoundCat infile outfile precision

   where
        infile:      input file
        outfile:     output file
        precision:   number of significant digits
   </pre>
*/
int main ( int argc, char ** argv)
{
  if ( argc != 4 ) usage(argc,argv);

  std::ifstream input(argv[1]);

  if ( !input )
  {
    std::cerr << "Error opening input file '" << argv[1] <<"'.\n";
    return 2;
  }

  input.unsetf(std::ios::skipws);

  std::ofstream output(argv[2]);

  if ( !output )
  {
    std::cerr << "Error opening output file '" << argv[2] <<"'.\n";
    return 3;
  }

  long int precision = std::strtol(argv[3], 0, 10); // returns zero upon failure
  if ( ! precision )
  {
    std::cerr << "Malformed precision argument: '" << argv[3] << "'\n";
    return 4;
  }

  output.precision(precision);

  while ( input )
  {
    double number;
    char letter;
    input >> number;
    bool not_a_number = input.fail();
    if ( not_a_number )
    {
      input.clear();
      input >> letter;
      output << letter;
    }
    else
    {
      output << number;
    }
  }
  return 0;
}
