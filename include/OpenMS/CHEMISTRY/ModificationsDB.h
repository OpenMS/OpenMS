// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_RESIDUEMODIFICATIONSDB_H
#define OPENMS_CHEMISTRY_RESIDUEMODIFICATIONSDB_H

#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CHEMISTRY/ResidueModification2.h>

#include <set>

namespace OpenMS
{
	// forward declarations
	class ResidueModification2;
	class Residue;

	/** @ingroup Chemistry
	
			@brief database which holds all residue modifications from unimod.org
			
			bla

			singleton
			
			...
	*/
	class ModificationsDB
	{					
		public:

			inline static ModificationsDB* getInstance()
      {
        static ModificationsDB* db_ = 0;
        if (db_ == 0)
        {
          db_ = new ModificationsDB;
        }
        return db_;
      }

			/// returns the number of modifications read from the unimod.xml file
			UInt getNumberOfModifications() const;

			/// 
			const ResidueModification2& getModification(UInt index) const throw(Exception::IndexOverflow);

			///
			
			
		protected:

			std::vector<ResidueModification2> mods_;
	
		private:

			/** @name Constructors and Destructors
      */
      //@{
      /// default constructor
      ModificationsDB();

      ///copy constructor
      ModificationsDB(const ModificationsDB& residue_db);

      /// constructor with filename where the residues are stored in
      ModificationsDB(const String& res_filename, const String& mod_filename)
        throw(Exception::FileNotFound, Exception::ParseError);

      /// destructor
      virtual ~ModificationsDB();
      //@}

      /** @name Assignment
      */
      //@{
      /// assignment operator
      ModificationsDB& operator = (const ModificationsDB& aa);
      //@}		
	};
}
#endif
