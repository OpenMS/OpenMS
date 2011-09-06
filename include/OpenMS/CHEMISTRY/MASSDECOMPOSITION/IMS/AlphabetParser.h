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
// $Maintainer: Stephan Aiche $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_ALPHABETPARSER_H
#define OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_ALPHABETPARSER_H

#include <fstream>
#include <istream>
#include <map>
#include <string>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS {

namespace ims {

/**
 * @brief An abstract templatized parser to load the data that is used to initialize @c Alphabet objects.
 *
 * @c AlphabetParser reads the input source, which is given as a template parameter @c InputSource , by
 * @c load (const std::string& fname) function where @c fname is the source name.
 * Loaded data can be retrieved by calling @c getElements().
 *
 * @see Alphabet
 */
template <typename AlphabetElementType = double, 
          typename Container = std::map<std::string, AlphabetElementType>,
          typename InputSource = std::istream>
class AlphabetParser
{
public:
  /**
   * Type of data to be loaded.
   */
  typedef Container ContainerType;

  /**
   * Loads the data from the InputSource with the name @c fname.
   * If there is an error occurred while reading data from InputSource,
   * @c IOException is thrown.
   *
   * @param fname The name of the input source.
   */
  void load(const std::string& fname);

  /**
   * Gets the data that was loaded.
   *
   * @return The data.
   */
  virtual ContainerType& getElements() = 0;

  /**
   * Parses the the given input source \c is \c.
   */
  virtual void parse(InputSource& is) = 0;

  /**
   * Destructor.
   */
  virtual ~AlphabetParser() {}
};

template <typename AlphabetElementType, typename Container, typename InputSource>
void AlphabetParser<AlphabetElementType, Container, InputSource>::load(const std::string& fname)
{
  std::ifstream ifs(fname.c_str());
  if (!ifs)
  {
    throw Exception::IOException(__FILE__, __LINE__, __PRETTY_FUNCTION__,fname);
  }
  this->parse(ifs);
}

} // namespace ims

} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_ALPHABETPARSER_H
