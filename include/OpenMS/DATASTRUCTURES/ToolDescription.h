// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Chris Bielow, Mathias Walzer $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_TOOLDESCRIPTION_H
#define OPENMS_DATASTRUCTURES_TOOLDESCRIPTION_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <map>

namespace OpenMS
{

	namespace Internal
	{

    struct FileMapping
    {
      String location; // a regex/macro mix; to be expanded by tool;
      String target;   // TOPP parameter that determines the desired name
       // thus: move location -> target

      FileMapping& operator=(const FileMapping& rhs)
      {
        if (this==&rhs) return *this;
        location = rhs.location;
        target = rhs.target;
        return *this;
      }
    };

    struct MappingParam
    {
      std::map<Int, String> mapping;
      FileMapping post_move;
      
      MappingParam& operator=(const MappingParam& rhs)
      {
        if (this==&rhs) return *this;
        mapping = rhs.mapping;
        post_move = rhs.post_move;
        return *this;
      }
    };

	  /**	
		  @brief ToolDescription Class.
  	
		  This class represents a ToolDescription.
  		
		  @ingroup Datastructures
	  */
	  struct ToolDescriptionGeneric
	  {
      bool is_internal;
      String name;
      String category;
      StringList type;

      ToolDescriptionGeneric& operator=(const ToolDescriptionGeneric& rhs)
      {
        if (this==&rhs) return *this;
        
        is_internal = rhs.is_internal;
        name = rhs.name;
        category = rhs.category;
        type = rhs.type;
        return *this;
      }
	  };

    struct ToolDescription :
      ToolDescriptionGeneric
    {
      String commandline;
      String path; //< filename to external tool
      MappingParam tr_table;
      Param param;

      ToolDescription& operator=(const ToolDescription& rhs)
      {
        if (this==&rhs) return *this;
        
        ToolDescriptionGeneric::operator=(rhs);
        commandline = rhs.commandline;
        path = rhs.path;
        tr_table = rhs.tr_table;
        param = rhs.param;
        return *this;
      }

    };
  
  } // namespace Internal

} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_TOOLDESCRIPTION_H
