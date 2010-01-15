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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_SOFTWARE_H
#define OPENMS_METADATA_SOFTWARE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/CVTermList.h>

namespace OpenMS 
{
	/**
		@brief Description of the software used for processing
		
		@ingroup Metadata
	*/
  class OPENMS_DLLAPI Software
  	: public CVTermList
  {
    public:
    	/// Constructor
      Software();
      /// Copy constructor
      Software(const Software& source);
      /// Destructor
      virtual ~Software();
 			
 			/// Assignment operator
      Software& operator= (const Software& source);

      /// Equality operator
      bool operator== (const Software& rhs) const;
      /// Equality operator
      bool operator!= (const Software& rhs) const;
			
			/// returns the name of the software
      const String& getName() const;
      /// sets the name of the software
      void setName(const String& name);
			
			/// returns the software version
      const String& getVersion() const;
      /// sets the software version
      void setVersion(const String& version);

    protected:
      String name_;
      String version_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_SOFTWARE_H
