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

#include <OpenMS/METADATA/MetaInfoDescription.h>

using namespace std;

namespace OpenMS
{
	
	MetaInfoDescription::MetaInfoDescription():
		MetaInfoInterface(),
		comment_(),
		name_()
	{
	  
	}
	
	MetaInfoDescription::MetaInfoDescription(const MetaInfoDescription& source):
		MetaInfoInterface(source),
	  comment_(source.comment_),
	  name_(source.name_),
	  source_file_(source.source_file_)
	{
	  
	}
	
	MetaInfoDescription::~MetaInfoDescription()
	{
	  
	}
	
	MetaInfoDescription& MetaInfoDescription::operator = (const MetaInfoDescription& source)
	{
	  if (&source == this) return *this;
	  
	  MetaInfoInterface::operator=(source);
	  comment_ = source.comment_;
	  name_ = source.name_;
	  source_file_ = source.source_file_;
	  
	  return *this;
	}

  bool MetaInfoDescription::operator== (const MetaInfoDescription& rhs) const
  {
  	return 
		  comment_ == rhs.comment_ &&
		  name_ == rhs.name_ &&
		  source_file_ == rhs.source_file_ &&
  		MetaInfoInterface::operator==(rhs)
  		;
  }
	
	const String& MetaInfoDescription::getComment() const 
	{
	  return comment_; 
	}
	
	void MetaInfoDescription::setComment(const String& comment)
	{
	  comment_ = comment; 
	}
	
	const SourceFile& MetaInfoDescription::getSourceFile() const 
	{
	  return source_file_; 
	}
	
	SourceFile&  MetaInfoDescription::getSourceFile()
	{
	  return source_file_; 
	}
	
	void MetaInfoDescription::setSourceFile(const SourceFile& source_file)
	{
	  source_file_ = source_file; 
	}

	void MetaInfoDescription::setName(const String& name)
	{
	  name_ = name; 
	}
	
}


