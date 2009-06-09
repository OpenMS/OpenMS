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
// $Maintainer: Guillaume Belz
// $Authors: Guillaume Belz
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/FidHandler.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <limits.h>

using namespace std;

namespace OpenMS
{
	namespace Internal
	{
      FidHandler::FidHandler(const String& filename)
        : ifstream(filename.c_str())
      {
        index_ = 0;
        intensityMin_ = INT_MAX;
        intensityMax_ = -INT_MAX;
        seekg(0, ios::beg);
			}

      unsigned int FidHandler::getIndex()
      {
        return index_;
      }
      
      unsigned int FidHandler::getIntensity()
      {
        // intensity is coded in 32 bits little-endian integer format
        unsigned char c1 = get();
        unsigned char c2 = get();
        unsigned char c3 = get();
        unsigned char c4 = get();
             
        unsigned int value = c4;
        value <<= 8;
        value |= c3;
        value <<= 8;
        value |= c2;
        value <<= 8;
        value |= c1;
    
        if(value > intensityMax_) intensityMax_ = value;
        if(value < intensityMin_) intensityMin_ = value;
        
        index_++;
        
        return value;
      }
      
      unsigned int FidHandler::getIntensityMin()
      {
        return intensityMin_;
      }

      unsigned int FidHandler::getIntensityMax()
      {
        return intensityMax_;
      }
      	
	} // namespace Internal
} // namespace OpenMS
