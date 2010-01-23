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

#include <OpenMS/FORMAT/DB/PersistentObject.h>

namespace OpenMS
{

	PersistentObject::PersistentObject():persistence_id_(0)
	{
		
	}

	PersistentObject::~PersistentObject()
	{
			
	}

	PersistentObject& PersistentObject::operator=(const PersistentObject& rhs)
	{
		if (&rhs==this) return *this;
		
		persistence_id_ = rhs.persistence_id_;
		
		return *this;
	}
	
  const UID& PersistentObject::getPersistenceId() const
	{
		return persistence_id_;
	}
	
  void PersistentObject::setPersistenceId(const UID& persistence_id)
	{
		persistence_id_ = persistence_id;
	}

	void PersistentObject::clearId(bool deep)
  {
    persistence_id_ = 0;
    if (deep == true)
    {
      clearChildIds_();
    }
  }

}
