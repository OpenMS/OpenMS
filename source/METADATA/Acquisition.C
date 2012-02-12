// -*- mode: C++; tab-width:s 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/METADATA/Acquisition.h>

using namespace std;

namespace OpenMS
{
	Acquisition::Acquisition():
		MetaInfoInterface(),
		identifier_()
	{
	}
	
	Acquisition::Acquisition(const Acquisition& source):
		MetaInfoInterface(source),
	  identifier_(source.identifier_)
	{
	}
	
	Acquisition::~Acquisition()
	{
	}
	
	Acquisition& Acquisition::operator = (const Acquisition& source)
	{
	  if (&source == this) return *this;
	  
	  identifier_ = source.identifier_;
	  MetaInfoInterface::operator=(source);
	  
	  return *this;
	}

  bool Acquisition::operator== (const Acquisition& rhs) const
  {
  	return 
  		identifier_ == rhs.identifier_ &&
  		MetaInfoInterface::operator==(rhs)
  		;
  }
  
  bool Acquisition::operator!= (const Acquisition& rhs) const
  {
  	return !(operator==(rhs));
 	}
	
	const String& Acquisition::getIdentifier() const 
	{
	  return identifier_; 
	}
	
	void Acquisition::setIdentifier(const String& identifier)
	{
	  identifier_ = identifier; 
	}

}	
	
