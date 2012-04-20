// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_ENHANCEDTABBARWIDGETINTERFACE_H
#define OPENMS_VISUAL_ENHANCEDTABBARWIDGETINTERFACE_H

#include <OpenMS/KERNEL/StandardTypes.h>

namespace OpenMS
{
  /**
    @brief Widgets that are placed into an EnhancedTabBar must implement this interface

    @ingroup Visual
  */
  class EnhancedTabBarWidgetInterface
  {
public:
    /// Destructor
    virtual ~EnhancedTabBarWidgetInterface() {}

    /// get the EnhancedTabBar unique window id
    virtual Int getWindowId() = 0;

    /// set the EnhancedTabBar unique window id
    virtual void setWindowId(Int window_id) = 0;
  };
}  // namespace OpenMS

#endif // OPENMS_VISUAL_ENHANCEDTABBARWIDGETINTERFACE_H
