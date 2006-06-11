// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: Software.h,v 1.1 2006/05/30 15:46:38 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_SOFTWARE_H
#define OPENMS_METADATA_SOFTWARE_H

#include <OpenMS/DATASTRUCTURES/String.h>

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
			
			/**
				@brief returns the duration of time used for processing
				
				The unit is minutes (default is 0).
			*/
      float getCompletionTime() const;
      /// sets the duration of time used for processing (in minutes)
      void setCompletionTime(float completion_time);


    protected:
      String name_;
      String version_;
      String comment_;
      float completion_time_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_SOFTWARE_H
