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
//
// --------------------------------------------------------------------------
// $Id: ModificationStringParser.h,v 1.4 2006/03/28 08:03:27 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------


#ifndef __MODIFICATIONSTRINGPARSER_H__
#define __MODIFICATIONSTRINGPARSER_H__

#include <string>
#include <ctype.h>                         //! isdigit(), isalpha()
#include <sstream>
#include <vector>


#include "Modification.h"


namespace OpenMS
  {

  /** used from class \c Sample to parse partial modifications out of strings:
   *  notation: position1 ( modification1 , modifcation2 , ... ) ; position2 ....
   *  the different modifications signify different possibilities for the position in question
   */
  class ModificationStringParser
    {
    private:

      //! reads a position out of \c temp_string
      bool readPosition(std::istringstream& temp_string_stream, int& position);

      //! reads modification out of \c temp_string
      bool readModification(std::istringstream& temp_string_stream, int& modification);

      //! database connection infos
      std::string db_username_, db_password_, db_host_;


    public:

      //! exception thrown if invalid formula string is given
    class InvalidModificationString : public OpenMS::Exception::Base
        {
        public:
          InvalidModificationString(const char* file, int line, const char* function, std::string mod_string) throw();
          ~InvalidModificationString() throw();
        };

      //! default constructor
      ModificationStringParser(std::string db_username, std::string db_password, std::string db_host);

      //! does the work
      std::vector<std::pair<int, std::vector<Modification*> > > parse(const std::string& mod_string) throw(InvalidModificationString);


    };
}

#endif



