// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/ToolDescription.h>

#include <map>

namespace OpenMS
{
  /**
     @brief Updates an INI
  */

  /// map each old TOPP tool to its new Name
  typedef std::map<Internal::ToolDescriptionInternal, Internal::ToolDescriptionInternal> ToolMapping;

  class OPENMS_DLLAPI INIUpdater
  {
public:

    INIUpdater();


    StringList getToolNamesFromINI(const Param & ini) const;

    const ToolMapping & getNameMapping();

    /*
      Finds the name of the new tool.
      The tools_type is optional and should be "" if there is none.

      The tools_type is ignored if there is a mapping without a type.

      @return true on success
    */
    bool getNewToolName(const String & old_name, const String & tools_type, String & new_name);

private:
    static ToolMapping map_;

  };

} // namespace OpenMS

