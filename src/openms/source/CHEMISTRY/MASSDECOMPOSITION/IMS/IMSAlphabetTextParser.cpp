// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
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
