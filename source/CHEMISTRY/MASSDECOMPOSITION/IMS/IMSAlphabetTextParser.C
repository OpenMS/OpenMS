// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> ?? $
// --------------------------------------------------------------------------
//

#include <sstream>

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabetTextParser.h>

/**
  Parses the data from the stream @c is .
  While loading the following is ignored:
    - white space
    - lines containing only white space
    - lines starting with '#' (even after leading whitespace, but not after anything else)

  @param is The input stream to be parsed.
*/
void OpenMS::ims::IMSAlphabetTextParser::parse(std::istream & is)
{
  // first make sure the store is empty
  elements_.clear();
  std::string line;
  std::string name;
  const std::string delimits(" \t"), comments("#");
  double mass;
  while (std::getline(is, line))
  {
    std::string::size_type i = line.find_first_not_of(delimits);
    if (i == std::string::npos || comments.find(line[i]) != std::string::npos)
    {
      continue;       // skip comment lines
    }
    std::istringstream input(line);
    input >> name >> mass;
    elements_.insert(std::make_pair(name, mass));
  }
}
