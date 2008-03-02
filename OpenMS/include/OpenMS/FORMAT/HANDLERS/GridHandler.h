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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_GRIDHANDLER_H
#define OPENMS_FORMAT_HANDLERS_GRIDHANDLER_H

#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <xercesc/sax2/Attributes.hpp>

namespace OpenMS
{
	class BaseMapping;
	class Grid;
	
  namespace Internal
  {
    /** 
    	@brief XML Handler for a vector of grid cells including their transformations.
    */
    class GridHandler
      : public XMLHandler
    {
	    public:
	      /// Constructor for loading 
	      GridHandler(Grid& grid, const String& filename);
	      /// Constructor dor storing
	      GridHandler(const Grid& grid, const String& filename);
	      /// Destructor
	      virtual ~GridHandler();
	
	      // Docu in base class
	      virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
	      // Docu in base class
	      virtual void characters(const XMLCh* const chars, unsigned int /*length*/);
	      // Docu in base class
	      virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
	      /// Print the contents to a stream
	      void writeTo(std::ostream& os);
	
	    protected:
	      /// Input grid
	      Grid* grid_;
	      /// Output grid
	      const Grid* cgrid_;
	      /// temporary variable for mapping type
	      BaseMapping* mapping_;
	      /// temporary parameters variable
	      Param param_;
				/// temporary variable for dimension
	      UInt dim_;

			private:
				///Not implemented
				GridHandler();
				
    }; // end of class GridHandler

  } // namespace Internal
} // namespace OpenMS

#endif
