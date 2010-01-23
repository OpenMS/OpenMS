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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_DB_PERSISTENTOBJECT_H
#define OPENMS_FORMAT_DB_PERSISTENTOBJECT_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <iostream>

namespace OpenMS
{ 
  /**
  	@brief Base class for all persistent objects.
  	
  	Interface for all classes that can be stored persistently in the OpenMS DB.

  	@ingroup DatabaseIO
  */
  class OPENMS_DLLAPI PersistentObject
  {

    public:
      /// Default constructor
      PersistentObject();

      /// Destructor
      virtual ~PersistentObject();
			
      /// Assignment operator
      PersistentObject& operator= (const PersistentObject& rhs);
			
      /**
      	@brief Returns the persistence id
      	
      	This id is only used in the DBAdapter the id is used to connect the object to the data stored in the DB.
      */
      const UID& getPersistenceId() const;

      /**
      	@brief Sets the persistence id
      	
      	This id is only used in the DBAdapter the id is used to connect the object to the data stored in the DB.
      	<BR>
      	Do not set the persistence id unless you know what you are doing!
      */
      void setPersistenceId(const UID& persistence_id);

      /**
      	@brief Clears the persistence id
      	
      	Sets the id to 0.<br>
      	@param deep determines which ids are cleared. <tt>false</tt> means that only the id of the current object is reset. 
      	<tt>true</tt> means that the ids of all sub-objects are reset as well (default).
      */
      void clearId(bool deep = true);

    protected:
    
      ///A persistence id used to refer the data back to the source
      UID persistence_id_;

      /**
      	@brief Clears the persistence id of all sub-objects.
      	
      	
      */
      virtual void clearChildIds_() =0;
  };

}
#endif
