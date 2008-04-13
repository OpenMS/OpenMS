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

#ifndef OPENMS_FORMAT_PEPNOVOINFILE_H
#define OPENMS_FORMAT_PEPNOVOINFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>

#include <map>


namespace OpenMS
{
	/**
		@brief PepNovo input file adapter.
		
		Creates a pepnovo.params file for PepNovo search from a peak list.
  	  	
  	@ingroup FileIO
	*/
  class PepNovoInfile
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

			/// stores the experiment data in a PepNovo input file that can be used as input for PepNovo shell execution
			String store(const String& filename) throw (Exception::UnableToCreateFile);

			/// retrieves the name, mass change, affected residues, type and position for all modifications from a string
			void handlePTMs(const String& modification_line, const String& modifications_filename, const bool monoisotopic) throw (Exception::FileNotReadable, Exception::FileNotFound, Exception::ParseError);

			/// return the modifications (the modification names map to the affected residues, the mass change and the type)
			const std::map< String, std::vector< String > >& getModifications() const;

		private:
			std::map< String, std::vector< String > > PTMname_residues_mass_type_;///< the modification names map to the affected residues, the mass change and the type
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_PEPNOVOINFILE_H
