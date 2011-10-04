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
#include <OpenMS/CONCEPT/LogStream.h>



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
      std::vector<FileMapping> pre_moves;
      std::vector<FileMapping> post_moves;
      
      MappingParam& operator=(const MappingParam& rhs)
      {
        if (this==&rhs) return *this;
        mapping = rhs.mapping;
        pre_moves = rhs.pre_moves;
        post_moves = rhs.post_moves;
        return *this;
      }
    };

	  /**	
		  @brief ToolDescription Class.
  	
		  This class represents a ToolDescription.
  		
		  @ingroup Datastructures
	  */
	  struct OPENMS_DLLAPI ToolDescriptionInternal
	  {
      bool is_internal;
      String name;
      String category;
      StringList types; // -types of the tool (if any, e.g. ['centroided','wavelet'])

      // default C'Tor
      ToolDescriptionInternal();

      // C'Tor with arguments
      ToolDescriptionInternal(const bool p_is_internal, const String& p_name, const String& p_category, const StringList& p_types);

      // short C'Tor
      ToolDescriptionInternal(const String& p_name, const StringList& p_types);

      ToolDescriptionInternal& operator=(const ToolDescriptionInternal& rhs);
      
      bool operator<(const ToolDescriptionInternal& rhs) const;
	  };

    struct OPENMS_DLLAPI ToolExternalDetails
    {
      String text_startup;
      String text_fail;
      String text_finish;
      String category;
      String commandline;
      String path; //< filename to external tool
      String working_directory; //< folder where the command will be executed from
      MappingParam tr_table;
      Param param;
    };

    /**
      Used for internal and external tools
    */
    struct OPENMS_DLLAPI ToolDescription :
      ToolDescriptionInternal
    {
      // additional details for external tools (one entry for each 'type')
      std::vector < ToolExternalDetails > external_details;

      // default CTor
      ToolDescription();

      // C'Tor for internal TOPP tools
      ToolDescription(const String& p_name, const String& p_category, const StringList& p_types = StringList());

      void addExternalType(const String& type, const ToolExternalDetails& details);

      void append(const ToolDescription& other);

      ToolDescription& operator=(const ToolDescription& rhs);
    };
  
  } // namespace Internal

} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_TOOLDESCRIPTION_H
