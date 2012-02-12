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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_PARAMXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_PARAMXMLHANDLER_H

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

#include <vector>
#include <map>

namespace OpenMS
{
	namespace Internal
	{
	/**
		@brief XML Handler for Param files.

	*/
  class OPENMS_DLLAPI ParamXMLHandler
  	: public XMLHandler
  {
    public:
    	/// Default constructor
      ParamXMLHandler(Param& param, const String& filename, const String& version);
			/// Destructor
      virtual ~ParamXMLHandler();

			// Docu in base class
      virtual void endElement( const XMLCh* const uri, const XMLCh* const local_name, const XMLCh* const qname);
			
			// Docu in base class
      virtual void startElement(const XMLCh* const uri, const XMLCh* const local_name, const XMLCh* const qname, const xercesc::Attributes& attributes);

    protected:
      /// The current absolute path (concatenation of nodes_ with <i>:</i> in between)
      String path_;
      /// Reference to the Param object to fill
      Param& param_;
      /// Map of node descriptions (they are set at the end of parsing)
			std::map<String,String> descriptions_;
      
      ///Temporary data for parsing of item lists
      struct
      {
        String name;
        String type;
        StringList stringlist;
        IntList intlist;
        DoubleList doublelist;
        StringList tags;
        String description;
        String restrictions;
        Int restrictions_index;
      } list_;

		private:
    	/// Not implemented
    	ParamXMLHandler();
  };

	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_PARAMXMLHANDLER_H
