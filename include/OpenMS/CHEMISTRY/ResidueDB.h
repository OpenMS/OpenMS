// -*- mode: C++; tab-width: 2; -*-
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

#ifndef OPENMS_CHEMISTRY_RESIDUEDB_H
#define OPENMS_CHEMISTRY_RESIDUEDB_H

#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <set>

namespace OpenMS
{
	// forward declarations
	class ResidueModification;
	class Residue;

	/** @ingroup Chemistry
	
			@brief residue data base which holds residues
			
			The residues stored in this DB are defined in a
			XML file under data/CHEMISTRY/residues.xml

			By default no modified residues are stored in an instance. However, if one
			queries the instance with getModifiedResidue, a new modified residue is 
			added. 
	*/
	class ResidueDB
	{					
		public:

			/** @name Typedefs
			*/
			//@{
			typedef std::set<Residue*>::iterator ResidueIterator;
			typedef std::set<const Residue*>::const_iterator ResidueConstIterator;
			//@}
			
			/// this member function serves as a replacement of the constructor
			inline static ResidueDB* getInstance()
      {
        static ResidueDB* db_ = 0;
        if (db_ == 0)
        {
          db_ = new ResidueDB;
        }
        return db_;
      }

			
			/** @name Constructors and Destructors
			*/
			//@{
			/// destructor
			virtual ~ResidueDB();
			//@}
			
			/** @name Accessors
			*/
			//@{
			/// returns the number of residues stored
			UInt getNumberOfResidues() const;

			UInt getNumberOfModifiedResidues() const;
			
			/// returns a pointer to the residue with name, 3 letter code or 1 letter code name
			const Residue* getResidue(const String& name) const;

			///
			const Residue* getModifiedResidue(const String& name);

			/// 
			const Residue* getModifiedResidue(const Residue* residue, const String& name);
			
			//const Residue* getModifiedResidue(const ResidueModification& mod);
			
			/// returns a set of all residues stored in this residue db
			const std::set<const Residue*>& getResidues() const;

			/// sets the residues from given file
			void setResidues(const String& filename);

			/// adds a residue, i.e. a unkown residue, where only the weight is known
			void addResidue(const Residue& residue);
			//@}

			/** @name Predicates
			*/
			//@{
			/// returns true if the db contains a residue with the given name
			bool hasResidue(const String& name) const;
			//@}
			
			/** @name Iterators
			*/
			//@{
			inline ResidueIterator beginResidue() { return residues_.begin(); }

			inline ResidueIterator endResidue() { return residues_.end(); }
			
			inline ResidueConstIterator beginResidue() const { return const_residues_.begin(); }

			inline ResidueConstIterator endResidue() const { return const_residues_.end(); }
			//@}

		protected:
      
			/** @name Private Constructors
			*/
			//@{
			/// default constructor
      ResidueDB();

      ///copy constructor
      ResidueDB(const ResidueDB& residue_db);
			//@}

			/** @name Assignment
      */
      //@{
      /// assignment operator
      ResidueDB& operator = (const ResidueDB& aa);
      //@}

			/// reads residues from the given file
			void readResiduesFromFile_(const String& filename);

			/// parses a residue, given the key/value pairs from i.e. an XML file
			Residue* parseResidue_(Map<String, String>& values) ;

			/// deletes all sub-instances of the stored data like modifications and residues
			void clear_();

			/// clears the residues
			void clearResidues_();

			/// builds an index of residue names for fast access, synonyms are also considered
			void buildResidueNames_();

			void addResidue_(Residue* residue);
			
			Map<String, Residue*> residue_names_;

			Map<String, Map<String, Residue*> > residue_mod_names_;

			std::set<Residue*> residues_;

			std::set<const Residue*> const_residues_;

			std::set<Residue*> modified_residues_;

			std::set<const Residue*> const_modified_residues_;
	};
}
#endif
