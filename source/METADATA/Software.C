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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Software.h>

using namespace std;

namespace OpenMS
{

	Software::Software()
		: CVTermList(),
	  	name_(),
	  	version_()
	{
	}
	
	Software::Software(const Software& rhs)
		: CVTermList(rhs),
	  	name_(rhs.name_),
	  	version_(rhs.version_)
	{
	}
	
	Software::~Software()
	{
	}
	
	Software& Software::operator = (const Software& rhs)
	{
	  if (&rhs == this) return *this;
	  
		CVTermList::operator = (rhs);
	  name_ = rhs.name_;
	  version_ = rhs.version_;
	  
	  return *this;
	}
	
	bool Software::operator== (const Software& rhs) const
	{
		return
			CVTermList::operator == (rhs) &&
	    name_ == rhs.name_ &&
	    version_ == rhs.version_
	    ;
	}
	
	bool Software::operator!= (const Software& rhs) const
	{
		return !(operator==(rhs));
	}
	 	
	const String& Software::getName() const 
	{
	  return name_; 
	}
	
	void Software::setName(const String& name)
	{
	  name_ = name; 
	}
	
	const String& Software::getVersion() const 
	{
	  return version_; 
	}
	
	void Software::setVersion(const String& version)
	{
	  version_ = version; 
	}

}

