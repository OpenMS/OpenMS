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
// $Id: Enzyme.h,v 1.4 2006/03/28 08:03:27 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------


#ifndef __ENZYME_H__
#define __ENZYME_H__


//STL-includes
#include <vector>
#include <list>
#include <string>
#include <ext/hash_map>

//OpenMS includes
#include <OpenMS/CONCEPT/Exception.h>

//MY includes
#include "string_hash_stl_fixes.h"
#include "../config_specannotate.h"
#include "AminoAcid.h"
#include "MySQLAdapter.h"


namespace OpenMS
  {

  //forward declarations
  class ProtDigMembers;


  //! This class represents enzymes used for cleaving the proteins
  class Enzyme
    {
    private:

      //! database information
      std::string db_username_, db_password_, db_host_;

      //! provides connection to MySQL database
      MySQLAdapter* sql_adapter_;

      //! unique ID of this enzyme in MySQL database
      int ID;

      //! gets unique ID for Enzyme out of Database, used in detailed constructor
      void getID_(std::string type) throw();


      //! contains type (name) of enzyme
      std::string type;

      //! contains cleavage sites
      std::vector<std::string> cleavage_sites;

      //! contains cleavage mode (C-terminal or N-terminal)
      std::string cleavage_mode;

    public:

      //! exception thrown if unknown enzyme is to be instanciated
    class UnknownEnzyme : public OpenMS::Exception::Base
        {
        public:
          UnknownEnzyme(const char* file, int line, const char* function, std::string request) throw();
          ~UnknownEnzyme() throw();
        };

      //! exception thrown if unknown cleavage mode is requested
    class UnknownCleavageMode : public OpenMS::Exception::Base
        {
        public:
          UnknownCleavageMode(const char* file, int line, const char* function, std::string request) throw();
          ~UnknownCleavageMode() throw();
        };

      //! constructor
      Enzyme();

      //! constructor with arguments: \c type specifies what enzyme is used in actual sample
      Enzyme(const std::string& ty, std::string db_username, std::string db_password,
             std::string db_host);

      //! destructor
      ~Enzyme();

      //! copy constructor
      Enzyme(const Enzyme& enzyme);

      /** gets necessary pointers to members of \c ProteinDigest, saves fragment start and stop indices in \c ProteinDigest::fragments
       * saves to each position all indices of \c ProteinDigest::fragments that signify fragments containing actual position in 
       * \c ProteinDigest::sequence_fragments
       */
      void digest(ProtDigMembers members) throw(UnknownCleavageMode);

      //! returns unique database ID of Enzyme
      int getID()
      {
        return ID;
      };

      //! returns type (name) of this Enzyme
      std::string getType()
      {
        return type;
      };


    };
}
#endif

