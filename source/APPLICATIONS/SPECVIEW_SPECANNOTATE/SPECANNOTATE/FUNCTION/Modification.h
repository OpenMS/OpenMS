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
// $Maintainer:  $
// --------------------------------------------------------------------------


#ifndef __MODIFICATION_H__
#define __MODIFICATION_H__

//STL includes
#include <string>
#include <algorithm>

//OpenMS includes
#include <OpenMS/CONCEPT/Exception.h>

//My includes
#include "../config_specannotate.h"
#include "MySQLAdapter.h"



namespace OpenMS
  {

  class ProtDigMembers;



  //! \todo evtl: ON DEMAND DATABASEACCESS FOR GETMASSES(-> Sample::calculateAnnotations schneller), destructor



  /** class representing modifications to certain residues
   * Integer accession codes for masses / formulae
   * 0: net formula / mass ADDED to molecule by modification
   * 1: net formula / mass SUBTRACTED from molecule by modification
   */
  class Modification
    {
    private:

      //! database information
      std::string db_username_, db_password_, db_host_;

      //! unique ID of this aminoacid in MySQL database
      int ID;

      //! provides connection to MySQL database
      MySQLAdapter* sql_adapter_;

      /** if connection to database is established this method fetches the unique database ID for actual instance of \c AminoAcid
       *    either from its one-letter-, three-letter-code or english name. result ist stored in member \c ID
       */
      void getID_(std::string type) throw();

      //! reads member values out of database
      void initialize_();

      //! contains the type (name) of modification
      std::string type;

      //! contains residues to which this modification should be applied
      std::vector<std::string> modification_sites;

      //! contains net molecular formula ADDED to molecule by modification
      std::string plus_formula;

      //! contains net molecular formula SUBTRACTED from molecule by modification
      std::string minus_formula;

      //! contains net monoisotopic molecular mass ADDED to molecule by modification
      std::string plus_mono_mass;

      //! contains net monoisotopic molecular mass SUBTRACTED from molecule by modification
      std::string minus_mono_mass;

      //! contains net average molecular mass ADDED to molecule by modification
      std::string plus_average_mass;

      //! contains net average molecular mass SUBTRACTED from molecule by modification
      std::string minus_average_mass;


    public:

      //! exception thrown if unknown position for aminoacid is requested
    class UnknownFormula : public OpenMS::Exception::Base
        {
        public:
          UnknownFormula(const char* file, int line, const char* function, int request) throw();
          ~UnknownFormula() throw();
        };

      //! exception thrown if unknown modification is to be instanciated
    class UnknownModification : public OpenMS::Exception::Base
        {
        public:
          UnknownModification(const char* file, int line, const char* function, std::string request) throw();
          ~UnknownModification() throw();
        };

      //! exception thrown if two overall modifications claim same residue
    class AmbiguousOverallModification : public OpenMS::Exception::Base
        {
        public:
          AmbiguousOverallModification(const char* file, int line, const char* function, std::string type) throw();
          ~AmbiguousOverallModification() throw();
        };

      //! constructor
      Modification();

      //! constructor with argument: \c type specifies the type of the modification
      Modification(const std::string& ty, std::string db_username, std::string db_password,
                   std::string db_host);

      //! constructor with argument: \c id specifies the database ID of the modification (if already known)
      Modification(int id, std::string db_username, std::string db_password,
                   std::string db_host);

      //! copy constructor
      Modification(const Modification& mod);

      //! modifies amino amino acids specified in \c modification_sites . gets pointers to necessarry members of \c ProteinDigest
      void modifyOverall(ProtDigMembers members);

      //! returns true, if this instance can modify given residue (given in one_letter_code)
      bool canModify(std::string residue);

      //! return type of instance of \c Modification
      std::string getType()
      {
        return type;
      };

      //! return database ID of instance of \c Modification
      int getID()
      {
        return ID;
      };

      //! this method changes ID (and type) of instance of modification, updates all member variables, WITHOUT establishing new datab. connect.
      void changeID(int new_ID);

      //! returns monoisotopic molecular mass \c sign: 0: + formula, 1 - formula
      float getMonoMass(int sign) throw();

      //! returns average molecular mass \c sign: 0: + formula, 1 - formula
      float getAverageMass(int sign) throw();

      //! returns molecular formula \c sign: 0: + formula, 1 - formula
      std::string getMolecularFormula(int sign) throw();
    };
}
#endif




