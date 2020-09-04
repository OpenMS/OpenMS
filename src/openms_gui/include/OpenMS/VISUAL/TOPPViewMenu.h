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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/DATASTRUCTURES/FlagSet.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/VISUAL/EnhancedWorkspace.h>

#include <QObject>

#include <vector>

class QAction;
class QMenu;

namespace OpenMS
{
  class TOPPViewBase;
  class EnhancedWorkspace;
  class RecentFilesMenu;
  
  enum class TV_STATUS
  {
    HAS_CANVAS,
    HAS_LAYER,
    HAS_MIRROR_MODE,  // implies 1D View
    IS_1D_VIEW,
    TOPP_IDLE
  };
  
  using FS_TV = FlagSet<TV_STATUS>;

  /// allow + operations on the enum, e.g. 'HAS_CANVAS + HAS_LAYER + IS_1D_VIEW'
  FS_TV OPENMS_GUI_DLLAPI operator+(const TV_STATUS left, const TV_STATUS right);

  /**
    @brief The file menu items for TOPPView



  */
  class TOPPViewMenu
    : public QObject
  {
    Q_OBJECT
  public:
    /** @brief Constructor which connects slots/signals of this class with the objects given as arguments

    @param parent Base class which actually shows the menu (as part of a QMainWindow)
    @param ws Workspace to connect some signals to
    @param recent_files A submenu for recent files which will be integrated as part of 'File -> Recent files'
    **/
    TOPPViewMenu(TOPPViewBase* const parent, EnhancedWorkspace* const ws, RecentFilesMenu* const recent_files);

    /// add a menu entry at 'Windows -> [Windowname]' to allow hiding/showing a TOPPView subwindow (e.g. Log, Layers, Filters, ...)
    void addWindowToggle(QAction* const window_toggle);

  public slots:
    /// enable/disable entries according to a given state of TOPPViewBase
    void update(const FS_TV status);

  private:
    struct ActionRequirement_
    {
      ActionRequirement_(QAction* action, const FS_TV& needs)
        : action(action), needs(needs) {}
      ActionRequirement_(QAction* action, const TV_STATUS& needs)
        : action(action), needs(needs) {}

      QAction* action;
      FS_TV needs;
    };

    /// fills menu_items_ members with ActionRequirements and returns the just created object
    /// Only use this for items which depend on the state of TOPPViewBase, 
    /// e.g. close() can only work if something is open. But open() is always allowed.
    const ActionRequirement_& addAction_(QAction* action, const TV_STATUS req);
    /// overload for multiple requirements
    const ActionRequirement_& addAction_(QAction* action, const FS_TV req);

    /// holds all actions which have a set of requirements, i.e. depend on the state of TOPPViewBase
    std::vector<ActionRequirement_> menu_items_;

    /// the windows submenu (holds all windows added via addWindowToggle())
    QMenu* m_windows_;
  };

} //namespace

