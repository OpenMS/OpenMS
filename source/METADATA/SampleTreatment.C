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

#include <OpenMS/METADATA/SampleTreatment.h>

using namespace std;

namespace OpenMS
{
		
	SampleTreatment::SampleTreatment():
		MetaInfoInterface(),
		comment_()
	{
		
	}
	
	SampleTreatment::SampleTreatment(const String& type):
		MetaInfoInterface(),
		type_(type),
		comment_()
	{
		
	}
	
	SampleTreatment::SampleTreatment(const SampleTreatment& source):
		MetaInfoInterface(source),
	  type_(source.type_),
		comment_(source.comment_)
	{
	  
	}
	
	SampleTreatment::~SampleTreatment()
	{
	  
	}
	
	SampleTreatment& SampleTreatment::operator = (const SampleTreatment& source)
	{
	  if (&source == this) return *this;
	  
	  MetaInfoInterface::operator=(source);
	  comment_ = source.comment_;
	  
	  return *this;
	}
	
	const String& SampleTreatment::getType() const
	{
		return type_;
	}
	
	const String& SampleTreatment::getComment() const
	{
		return comment_;
	}
	
	void SampleTreatment::setComment(const String& comment)
	{
		comment_ = comment;
	}
	
	bool SampleTreatment::operator== (const SampleTreatment& rhs) const
	{
		return
			MetaInfoInterface::operator==(rhs) &&
			comment_ == rhs.comment_
			;
	}
}

