// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Nico Pfeifer $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PEPXMLFILEMASCOT_H
#define OPENMS_FORMAT_PEPXMLFILEMASCOT_H

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

#include <vector>

namespace OpenMS 
{
  /**
    @brief Used to load Mascot PepXML files
    
		A schema for this format can be found at http://www.matrixscience.com/xmlns/schema/pepXML_v18/pepXML_v18.xsd. 
		  	
  	@ingroup FileIO
  */
  class OPENMS_DLLAPI PepXMLFileMascot
  	: protected Internal::XMLHandler,
  		public Internal::XMLFile
  {
		public:
		
			/// Constructor
			PepXMLFileMascot();
			
			/**
				@brief Loads peptide sequences with modifications out of a PepXML file

				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
			*/
			void load(const String& filename,  std::map<String, std::vector<AASequence> >& peptides);

  	protected:
  		
			// Docu in base class
			virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
			
			// Docu in base class
			virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
			
		  void matchModification_(DoubleReal mass, String& modification_description);
		  												
			/// @name members for loading data
			//@{
			/// Pointer to fill in protein identifications

			/// The title of the actual spectrum
			String actual_title_;
			
			/// The sequence of the actual peptide hit				
			String actual_sequence_;
			
			/// The modifications of the actual peptide hit (position is 1-based)
			std::vector<std::pair<String, UInt> > actual_modifications_;
			
			/// The peptides together with the spectrum title
			std::map<String, std::vector<AASequence> >* peptides_;
			
			/// stores the actual peptide sequences
			std::vector<AASequence> actual_aa_sequences_;
			
			/// stores the fixed residue modifications 
			std::vector<String> fixed_modifications_;

			/// stores the variable residue modifications
			std::vector<std::pair<String, DoubleReal> > variable_modifications_;
			//@}
  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_PEPXMLFILEMASCOT_H
