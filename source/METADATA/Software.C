// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/METADATA/Software.h>

using namespace std;

namespace OpenMS
{

	Software::Software():
	  name_(),
	  version_(),
	  comment_(),
		completion_time_()
	{
	  
	}
	
	Software::Software(const Software& source):
	  name_(source.name_),
	  version_(source.version_),
	  comment_(source.comment_),
	  completion_time_(source.completion_time_)
	{
	  
	}
	
	Software::~Software()
	{
	  
	}
	
	Software& Software::operator = (const Software& source)
	{
	  if (&source == this) return *this;
	  
	  name_ = source.name_;
	  version_ = source.version_;
	  comment_ = source.comment_;
	  completion_time_ = source.completion_time_;
	  
	  return *this;
	}
	
	bool Software::operator== (const Software& rhs) const
	{
		return 
	    name_ == rhs.name_ &&
	    version_ == rhs.version_ &&
	    comment_ == rhs.comment_ &&
	    completion_time_ == rhs.completion_time_
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
	
	const String& Software::getComment() const 
	{
	  return comment_; 
	}
	
	void Software::setComment(const String& comment)
	{
	  comment_ = comment; 
	}
	
	const DateTime& Software::getCompletionTime() const
	{
	  return completion_time_;
	}

	void Software::setCompletionTime(const DateTime& completion_time)
	{
	  completion_time_ = completion_time;
	}

	void Software::setCompletionTime(const String& completion_time)
	{
		completion_time_.set(completion_time);
	}

}

