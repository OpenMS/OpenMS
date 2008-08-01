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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_XMLVALIDATOR_H
#define OPENMS_FORMAT_XMLVALIDATOR_H

#include <OpenMS/DATASTRUCTURES/String.h>

#include <xercesc/sax/ErrorHandler.hpp>

namespace OpenMS
{
	/**
		@brief Validator for XML files.
		
		Validates an XML file against a given schema.
		
  	@ingroup FileIO
	*/
  class XMLValidator
  	: private xercesc::ErrorHandler
  {
    public:
    	/// Constructor
    	XMLValidator();

			/**
			  @brief Returns if an XML file is valid for given a schema file

			  @exception Exception::FileNotFound is thrown if the file cannot be found
			  @exception Exception::ParseError is thrown if the parser could not be initialized
			*/
			bool isValid(const String& filename, const String& schema);

  	protected:
  		//
  		bool valid_;
  		
  		/// @name Implementation of Xerces ErrorHandler methods
  		//@{
  		virtual void warning(const xercesc::SAXParseException& exception);
			virtual void error(const xercesc::SAXParseException& exception);
			virtual void fatalError(const xercesc::SAXParseException& exception);
			virtual void resetErrors();
  		//@}
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_XMLVALIDATOR_H
