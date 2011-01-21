// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_TRANSFORMATIONXMLFILE_H
#define OPENMS_FORMAT_TRANSFORMATIONXMLFILE_H


#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>

#include <vector>

namespace OpenMS 
{
  /**
	@brief Used to load and store TransformationXML files
    
	This class is used to load and store documents that implement the schema of
	TransformationXML files.
		
	A documented schema for this format can be found at http://open-ms.sourceforge.net/schemas/.
		
	@ingroup FileIO
  */
  class OPENMS_DLLAPI TransformationXMLFile
  	: protected Internal::XMLHandler,
  		public Internal::XMLFile
  {
	 public:
		
		/// Constructor
		TransformationXMLFile();
			
		/**
		@brief Loads the transformation from an TransformationXML file
				
		The information is read in and the information is stored in the
		corresponding variables

		@exception Exception::FileNotFound is thrown if the file could not be opened
		@exception Exception::ParseError is thrown if an error occurs during parsing
		*/
		void load(const String& filename, TransformationDescription& transformation);
			 			 
		/**
		@brief Stores the data in an TransformationXML file
				
		The data is read in and stored in the file named 'filename'.

		@exception Exception::UnableToCreateFile is thrown if the file could not be created
		@exception Exception::IllegalArgument is thrown if unsupported parameter types have been set
		*/
		void store(String filename, const TransformationDescription& transformation); 
  	
	 protected:
		// Docu in base class
		virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
			
		/// @name Members for use during loading data
		//@{
		/// Param to fill in during parse
		Param params_;
		/// Data vector
		TransformationDescription::DataPoints data_;
		/// Model type
		String model_type_;
		//@}

  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_TRANSFORMATIONXMLFILE_H
