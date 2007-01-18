// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/DATASTRUCTURES/HashMap.h>

#include <set>
#include <vector>

namespace OpenMS
{
	/** @ingroup Chemistry
	
			@brief residue data base which holds residues
			
			The residues stored in this DB are defined in a
			XML file under data/CHEMISTRY/residues.xml
	*/
	class ResidueDB
	{					
		public:

			/** @name Typedefs
			*/
			//@{
			typedef std::set<Residue*>::iterator ResidueIterator;
			//typedef std::set<const Residue*>::const_iterator ResidueConstIterator;
			typedef std::set<Residue::Modification*>::iterator ModificationIterator;
			//typedef std::set<const Residue::Modification*>::const_iterator ModificationConstIterator;
			//@}
			
			/** @name Constructors and Destructors
			*/
			//@{
			/// default constructor
			ResidueDB();

			///copy constructor
			ResidueDB(const ResidueDB& residue_db);
			
			/// constructor with filename where the residues are stored in
			ResidueDB(const String& res_filename, const String& mod_filename) 
				throw(Exception::FileNotFound, Exception::ParseError);
			
			/// destructor
			virtual ~ResidueDB();
			//@}
			
			/** @name Assignment
			*/
			//@{
			/// assignment operator
			ResidueDB& operator = (const ResidueDB& aa);
			//@}
			
			/** @name Accessors
			*/
			//@{
			/// returns the number of residues stored
			Size getNumberOfResidues() const;
			
			/// returns the number of modifications stored in this residue db
			Size getNumberOfModifications() const;

			/// resturns a pointer to modification with name name, if non is found 0 is returned
			const Residue::Modification* getModification(const String& name) const;

			/// returns a set of modifications which can be applied to the given residue
			std::set<const Residue::Modification*> getModifications(const Residue* residue) const;

			/// returns a set of modifications which can be applied to the given residue
			std::set<const Residue::Modification*> getModifications(const String& res_name) const;

			/// returns a set of all modifications stored in this residue db
			const std::set<const Residue::Modification*>& getModifications() const;

			/// returns a pointer to the residue with name, 3 letter code or 1 letter code name
			const Residue* getResidue(const String& name) const;

			/// returns a set of residues which can have the given modification 
			std::set<const Residue*> getResidues(const Residue::Modification* modification) const;

			/// returns a set of residues which can have the given modification
			std::set<const Residue*> getResidues(const String& mod_name) const;

			/// returns a set of all residues stored in this residue db
			const std::set<const Residue*>& getResidues() const;

			/// sets the modifications from given file
			void setModifications(const String& filename) throw(Exception::FileNotFound, Exception::ParseError);

			/// adds a modification, i.e. an unknown modification, where only the weights are known
			void addModification(Residue::Modification modification);

			/// sets the residues from given file
			void setResidues(const String& filename) throw(Exception::FileNotFound, Exception::ParseError);

			/// adds a residue, i.e. a unkown residue, where only the weight is known
			void addResidue(Residue residue);
			//@}

			/** @name Predicates
			*/
			//@{
			/// returns true if the db contains a modification with the given name
			bool hasModification(const String& name) const;

			/// returns true if the db contains a residue with the given name
			bool hasResidue(const String& name) const;
			//@}
			
			/** @name Iterators
			*/
			//@{
			inline ResidueIterator beginResidue() { return residues_.begin(); }

			inline ResidueIterator endResidue() { return residues_.end(); }
			
			//inline ResidueConstIterator beginResidue() const { return const_residues_.begin(); }

			//inline ResidueConstIterator endResidue() const { return const_residues_.end(); }

			inline ModificationIterator beginModification() { return modifications_.begin(); }

			inline ModificationIterator endModification() { return modifications_.end(); }
				
			//inline ModificationConstIterator beginModification() const { return const_modifications_.begin(); }

			//inline ModificationConstIterator endModification() const { return const_modifications_.end(); }
			//@}

		protected:

			/*_ reads residues from the given file
			*/
			void readResiduesFromFile_(const String& filename) throw(Exception::FileNotFound, Exception::ParseError);

			/*_ parses a residue, given the key/value pairs from i.e. an XML file
			*/
			Residue* parseResidue_(HashMap<String, String>& values) throw();

			/*_ reads modifications from a file
			*/
			void readModificationsFromFile_(const String& filename) throw(Exception::FileNotFound, Exception::ParseError);

			/*_ deletes all sub-instances of the stored data like modifications and residues
			*/
			void clear_();

			/*_ deletes all residues
			*/
			void clearResidues_();

			/*_ deletes all modifications and also modified residues
			*/
			void clearModifications_();

			/*_ builds an index of residue names for fast access, synonyms are also considered
			*/
			void buildResidueNames_();

			/*_ builds an index of modifications names for fast access, synonyms are also considered
			*/
			void buildModificationNames_();

			/*_ builds modified residues from a given modifications and residues
			*/
			void buildModifiedResidues_();

			HashMap<String, Residue*> residue_names_;

			std::set<Residue*> residues_;

			std::set<const Residue*> const_residues_;
		
			HashMap<String, Residue::Modification*> modification_names_;
			
			std::set<Residue::Modification*> modifications_;

			std::set<const Residue::Modification*> const_modifications_;
	};
}
#endif
