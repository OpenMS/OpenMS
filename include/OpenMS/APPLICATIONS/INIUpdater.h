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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_APPLICATIONS_INIUPDATER_H
#define OPENMS_APPLICATIONS_INIUPDATER_H

#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/ToolDescription.h>

namespace OpenMS
{
  /**
  	 @brief Updates an INI

     


  */
  
  /// map each old TOPP/UTIL to its new Name
  typedef Map<Internal::ToolDescriptionInternal, Internal::ToolDescriptionInternal> ToolMapping;

  class OPENMS_DLLAPI INIUpdater
  {
    public:
   
      INIUpdater();
      

      StringList getToolNamesFromINI(const Param& ini) const;
    
      const ToolMapping& getNameMapping();

      /*
        Finds the name of the new tool.
        The tools_type is optional and should be "" if there is none.

        The tools_type is ignored if there is a mapping without a type.

        @return true on success
      */
      bool getNewToolName(const String& old_name, const String& tools_type, String& new_name);

  private:
      static ToolMapping map_;

  };

} // namespace OpenMS

#endif //OPENMS_APPLICATIONS_INIUPDATER_H
