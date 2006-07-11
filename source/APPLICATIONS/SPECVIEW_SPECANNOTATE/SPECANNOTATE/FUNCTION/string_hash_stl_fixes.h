// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// --------------------------------------------------------------------------


#ifndef __STRING_HASH_STL_FIXES_H__
#define __STRING_HASH_STL_FIXES_H__


#include <string>
#include <OpenMS/DATASTRUCTURES/String.h>
//#include <ext/stl_hash_fun.h>

using namespace OpenMS;

namespace __gnu_cxx
  {
  //fix for hashing stl strings with sgi hash_map
  template<>
  struct hash<std::string>
    {
      size_t operator()(const std::string& x) const
        {
          return hash<const char*>()( x.c_str() );
        }
    };

  //fix for hashing Strings with sgi hash_map
  template<>
  struct hash<String>
    {
      size_t operator()(const String& x) const
        {
          return hash<const char*>()( x.c_str() );
        }
    };
}


#endif
