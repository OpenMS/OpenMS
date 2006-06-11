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
// $Id: AminoAcid.h,v 1.4 2006/03/28 08:03:27 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------


#ifndef __AMINOACID_H__
#define __AMINOACID_H__

//STL includes
#include <string>
#include <map>

//OpenMS includes
#include <OpenMS/CONCEPT/Exception.h>

//My includes
#include "../config_specannotate.h"
#include "MySQLAdapter.h"



namespace OpenMS
  {
  //! class representing amino acids.
  class AminoAcid
    {
    private:

      //! database information
      std::string db_username_, db_password_, db_host_;

      //! unique ID of this aminoacid in MySQL database
      int ID;

#ifdef ANNOTATE_XML
      //! maps one-letter-code and two_letter_code to names and also names to themselves
      std::map<std::string, std::string> code_names_;
#endif

      //! contains aminoacid name
      std::string name_;

      //! contains one-letter-code of aminoacid
      std::string one_letter_code_;

      //! contains three-letter-code of aminoacid
      std::string three_letter_code_;

      //! contains monoisotopic masses
      double middle_mono_mass_;
      double n_term_mono_mass_;
      double c_term_mono_mass_;
      double single_mono_mass_;

      //! contains average masses
      double middle_average_mass_;
      double n_term_average_mass_;
      double c_term_average_mass_;
      double single_average_mass_;

      //! contains molecular formulae
      std::string middle_formula_;
      std::string n_term_formula_;
      std::string c_term_formula_;
      std::string single_formula_;

      //! provides connection to MySQL database
      MySQLAdapter* sql_adapter_;

      /** if connection to database is established this method fetches the unique database ID for actual instance of \c AminoAcid
       *    either from its one-letter-, three-letter-code or english name. result ist stored in member \c ID
       */
      void getID_(std::string type) throw();

      //! this method connects to database calls getID_ and then fills member variables out of database
      void initialize_(std::string type);

    public:

      //! exception thrown if unknown position for aminoacid is requested
    class UnknownPosition : public OpenMS::Exception::Base
        {
        public:
          UnknownPosition(const char* file, int line, const char* function, int request) throw();
          ~UnknownPosition() throw();
        };

      //! exception thrown if unknown aminoacid is to be instanciated
    class UnknownAminoAcid : public OpenMS::Exception::Base
        {
        public:
          UnknownAminoAcid(const char* file, int line, const char* function, std::string request) throw();
          ~UnknownAminoAcid() throw();
        };


      //! default constructor
      AminoAcid();

      /** constructor that is initialized
       *    by a string containing type of amino acid in one- or three-letter-code, or the english name
       *    and by some strings containing data for database connection
       */
      AminoAcid(std::string type, std::string db_username, std::string db_password,
                std::string db_host);

      //! destructor
      ~AminoAcid();

      //! copy constructor
      AminoAcid(const AminoAcid& amino_acid);

      //! assignment operator
      const AminoAcid& operator=(const AminoAcid& amino_acid);

      //! returns the formula as string: @param position: 0: middle, 1: C-terminal, 2: N-terminal, 3: single
      std::string getFormula(int position) throw(UnknownPosition);

      //! returns the monoisotopic mass: @param position: 0: middle, 1: C-terminal, 2: N-terminal, 3: single
      double getMonoMass(int position) throw(UnknownPosition);

      //! returns the average mass: @param position: 0: middle, 1: C-terminal, 2: N-terminal, 3: single
      double getAverageMass(int position) throw(UnknownPosition);

      //! returns name of amino acid
      std::string getName() const;

      //! returns one-letter-code of amino acid
      std::string getOneLetter() const;

      //! returns three-letter-code of amino acid
      std::string getThreeLetter() const;

      //! returns the ID
      int getID()
      {
        return ID;
      };

    };

}
#endif //__AMINOACID_H__
