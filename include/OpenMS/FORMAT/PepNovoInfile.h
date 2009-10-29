// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
		
		Creates a pepnovo.params file for PepNovo search from a peak list.
  	  	
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

			/** retrieves the name, mass change, affected residues, type and position for all modifications from a string
					
					@param modification_line
					@param modifications_filename
					@param monoisotopic
					@throw FileNotFound is thrown if the given file is not found
					@throw FileNotReadable is thrown if the given file is not readable
					@throw ParseError is throw if the given file could not be parsed
			*/
			String handlePTMs(const String &modification, const bool variable);

			void setModifications(const StringList &fixed_mods, const StringList &variable_mods);

			/// return the modifications (the modification names map to the affected residues, the mass change and the type)
			///void getModifications(std::vector<Modification>&) const;
			void getModifications(std::map<String,String>& modification_key_map) const;

		private:
			ModificationDefinitionsSet mods_;
			std::map<String,String>mods_and_keys_;
			TextFile ptm_file_;
			void generate_pepnovo_lines();

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_PEPNOVOINFILE_H
