// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Sandro Andreotti $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PEPNOVOINFILE_H
#define OPENMS_FORMAT_PEPNOVOINFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <map>


namespace OpenMS
{
	/**
		@brief PepNovo input file adapter.
		
		Creates a PepNovo_PTMs.txt file for PepNovo search.
  	  	
  	@ingroup FileIO
	*/
  class OPENMS_DLLAPI PepNovoInfile
  {
		public:
			/// default constructor
			PepNovoInfile();

			/// copy constructor
			PepNovoInfile(const PepNovoInfile& pepnovo_infile);

			/// destructor
			virtual ~PepNovoInfile();

			/// assignment operator
			PepNovoInfile& operator=(const PepNovoInfile& pepnovo_infile);

			/// equality operator
			bool operator==(const PepNovoInfile& pepnovo_infile) const;

			/** stores the experiment data in a PepNovo input file that can be used as input for PepNovo shell execution
					
				@param filename the file which the input file is stored into
				@throw Exception::UnableToCreateFile is thrown if the given file could not be created
			*/
			void store(const String& filename);

			/** @brief generates the PepNovo Infile for given fixed and variable modifications			 *
			 *
			 * @param fixed_mods StringList of fixed modifications unique identifiers
			 * @param variable_mods StringList of variable modifications unique identifiers
			 */
			void setModifications(const StringList &fixed_mods, const StringList &variable_mods);

			/** @brief return the modifications.
			 *
			 *  the modification unique identifiers are mapped to the keys used
			 *  in the PepNovo Infile (origin+rounded monoisotopic mass of modification ).
			 *  (e.g. modification_key_map["K+16"]=="Oxidation (K)" )
			 */
			void getModifications(std::map<String,String>& modification_key_map) const;

		private:
			ModificationDefinitionsSet mods_;
			std::map<String,String>mods_and_keys_;
			TextFile ptm_file_;


		 /** retrieves the name of modification, and generates the corresponding line for the
			 PepNovo infile.
			 @param modification the modification
			 @param variable should be set to true if it variable
		*/
		  String handlePTMs_(const String &modification, const bool variable);
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_PEPNOVOINFILE_H
