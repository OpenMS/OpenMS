// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                    Flex series file support
// --------------------------------------------------------------------------
//  Copyright (C) 2009 -- Guillaume Belz (guillaume.belz@chu-lyon.fr)
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
// $Maintainer: Guillaume Belz$
// $Authors: Guillaume Belz$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/FidHandler.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <climits>

using namespace std;

namespace OpenMS
{
	namespace Internal
	{
      FidHandler::FidHandler(const String& filename)
        : ifstream(filename.c_str(), ios_base::binary | ios_base::in)
      {
        index_ = 0;
        seekg(0, ios::beg);
			}

      FidHandler::~FidHandler()
      {
      }
      
      Size FidHandler::getIndex()
      {
        return index_;
      }
      
      Size FidHandler::getIntensity()
      {
        // intensity is coded in 32 bits little-endian integer format
        Size c1 = get();
        Size c2 = get();
        Size c3 = get();
        Size c4 = get();

        Size value = c4;
        value <<= 8;
        value |= c3;
        value <<= 8;
        value |= c2;
        value <<= 8;
        value |= c1;
    
        index_++;
        
        return value;
      }
	} // namespace Internal
} // namespace OpenMS

