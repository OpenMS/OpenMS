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
#include <OpenMS/CHEMISTRY/ResidueModification.h>

#include <set>

namespace OpenMS
{
	// forward declarations
	class ResidueModification;
	class Residue;

	/** @ingroup Chemistry
	
			@brief database which holds all residue modifications from PSI-MOD
			
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

			/// returns the modification with the given index 
			const ResidueModification& getModification(UInt index) const;

			/** @brief returns the modification of the given name

					This can either be the PSI-MOD identifier or every other unique
					identifier which can be found in the PSI-MOD definitions file.

			*/
			const ResidueModification& getModification(const String& name) const;

			void getModificationsByDiffMonoMass(std::vector<String>& mods, double mass, double error = 0.0);
			
			void getModificationsByDiffMonoMass(std::vector<String>& mods, const String& residue, double mass, double error = 0.0);


			/// adds modifications from a given file in OBO format
			void readFromOBOFile(const String& filename);

			/// adds modifications from a given file in Unimod XML format
			void readFromUnimodXMLFile(const String& filename);
			
		protected:

			/// stores the modifications
			std::vector<ResidueModification*> mods_;

			/// stores the mappings of (unique) names to the modifications
			Map<String, const ResidueModification*> modification_names_;
			
			
		private:

			/** @name Constructors and Destructors
      */
      //@{
      /// default constructor
      ModificationsDB();

      ///copy constructor
      ModificationsDB(const ModificationsDB& residue_db);

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
