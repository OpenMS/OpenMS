// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

namespace OpenMS
{
  class EnhancedTabBar;
  /**
    @brief Widgets that are placed into an EnhancedTabBar must implement this interface

    @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI EnhancedTabBarWidgetInterface
  {
public:
    /// C'tor; creates a new ID;
    EnhancedTabBarWidgetInterface();
    /// Destructor
    virtual ~EnhancedTabBarWidgetInterface();

    /// adds itself to this tabbar and upon destruction removes itself again.
    /// Make sure the tabbar still exists when you call this function and this object is destroyed
    void addToTabBar(EnhancedTabBar* const parent, const String& caption, const bool make_active_tab = true);

    /// get the EnhancedTabBar unique window id
    Int getWindowId();

    /// the first object to be created will get this ID
    static Int getFirstWindowID();

private:
    
    Int window_id_ { -1 };
    EnhancedTabBar* parent_ { nullptr }; ///< this is our parent. Which will delete us automatically when we destruct
  };
}  // namespace OpenMS

