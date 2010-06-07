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

#ifdef OPENMS_BIG_ENDIAN
  template<typename T> T ByteReverse(const T in)
  {
    T out;
    const char* pin = (const char*) &in;
    char* pout= (char*) (&out+1) - 1 ;

    int i;
    for(i= sizeof(T) ; i>0 ; --i)
    {
      *pout-- = *pin++ ;
    }
    return out ;
  }
#endif

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
        Int32 result = 0;
        read( (char*) &result, 4);
        #ifdef OPENMS_BIG_ENDIAN
          result = ByteReverse<Int32>(result);
        #endif
        index_++;        
        return (result > 0) ? result : 0;
      }
	} // namespace Internal
} // namespace OpenMS
