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

#ifndef OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_IMSALPHABETTEXTPARSER_H
#define OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_IMSALPHABETTEXTPARSER_H

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabetParser.h>

namespace OpenMS {

  namespace ims {

    /**
      @brief Implements abstract @c AlphabetParser to read data from the plain text format.

      @c AlphabetTextParser parses the data source using overriden @c parse(std::istream&)
      and stores the parsed data permanently. That can be retrieved by @c getElements() function.
    */
    class OPENMS_DLLAPI IMSAlphabetTextParser 
      : public IMSAlphabetParser<> 
    {
    private:
      /**
        The parsed data.
      */
      ContainerType elements_;
    public:
      /**
        Gets the parsed data.

        @return The parsed data.
      */
      virtual ContainerType& getElements() { return elements_; }

      /**
        Parses the input stream \c is \c.

        @param is The input stream to be parsed
      */
      virtual void parse(std::istream& is);
    };

  } // namespace ims
} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_ALPHABETTEXTPARSER_H
