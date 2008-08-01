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

#ifndef OPENMS_METADATA_SOFTWARE_H
#define OPENMS_METADATA_SOFTWARE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

namespace OpenMS 
{
	/**
		@brief Description of the software used for processing
		
		 
		
		@ingroup Metadata
	*/
  class Software
  {
    public:
    	/// Constructor
      Software();
      /// Copy constructor
      Software(const Software& source);
      /// Destructor
      ~Software();
 			
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
			
			/// returns the free-text comment
      const String& getComment() const;
      /// sets the free-text comment
      void setComment(const String& comment);

			/// returns the time of completition of the processing
	    const DateTime& getCompletionTime() const;
      /// sets the time of completition taking a DateTime object
      void setCompletionTime(const DateTime& completion_time);
      /// sets the time of completition taking a String object
			/// provided for convenience
      void setCompletionTime(const String& completion_time);


    protected:
      String name_;
      String version_;
      String comment_;
      DateTime completion_time_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_SOFTWARE_H
