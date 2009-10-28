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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/UniqueIdInterface.h>
#include <iomanip>

namespace OpenMS
{

void
UniqueIdInterface::setUniqueId(const String & rhs)
{
  clearUniqueId();

  String::size_type last_underscore = rhs.rfind('_');
  // Note: String::npos is usually defined as size_t(-1).
  // In any case, the condition of the next if() statement evaluates to a constant
  // and can be "compiled away".
  if ( String::npos + 1 != 0 )
  {
    if ( last_underscore == String::npos )
    {
      last_underscore = -1;
    }
  }
  // For the next line to be correct in case rhs contains no '_', it is necessary that npos+1==0;
	String s=rhs.substr(last_underscore + 1);

	for (String::const_iterator s_i=s.begin();s_i<s.end();++s_i)
	{
		int i = (*s_i -'0');
		if (i<0 || i>9) 
		{
			clearUniqueId();
			return;
		}
    unique_id_ = 10*unique_id_ + i;
  }
	
}

}
