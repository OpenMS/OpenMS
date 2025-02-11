// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
#include <OpenMS/VISUAL/LayerDataBase.h>

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

  using FS_LAYER = FlagSet<LayerDataBase::DataType>;
  /// allow + operations on the enum, e.g. 'DT_PEAK + DT_FEATURE'
  FS_LAYER OPENMS_GUI_DLLAPI operator+(const LayerDataBase::DataType left, const LayerDataBase::DataType right);


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
    void update(const FS_TV status, const LayerDataBase::DataType layer_type);

  private:
    struct ActionRequirement_
    {
      ActionRequirement_(QAction* action, const FS_TV& needs, const FS_LAYER layer_set)
        : action_(action), needs_(needs), layer_set_(layer_set) {}
      ActionRequirement_(QAction* action, const TV_STATUS& needs, const FS_LAYER layer_set)
        : action_(action), needs_(needs), layer_set_(layer_set) {}

      /// check if an ActionRequirement is fulfilled by the arguments
      /// i.e. @p status is a superset of needs_ and @p layer_type is a superset of layer_set_ (or layer_set_ is empty)
      /// If yes, the the action to enabled, or disable it otherwise
      void enableAction(const FS_TV status, const LayerDataBase::DataType layer_type);

    private:
      QAction* action_;
      FS_TV needs_;
      FS_LAYER layer_set_;
    };

    /// fills menu_items_ members with ActionRequirements and returns the just created object
    /// Only use this for items which depend on the state of TOPPViewBase, 
    /// e.g. close() can only work if something is open. But open() is always allowed.
    QAction* addAction_(QAction* action, const TV_STATUS req, const FS_LAYER layer_set = FS_LAYER());
    /// overload for multiple requirements
    QAction* addAction_(QAction* action, const FS_TV req, const FS_LAYER layer_set = FS_LAYER());

    /// holds all actions which have a set of requirements, i.e. depend on the state of TOPPViewBase
    std::vector<ActionRequirement_> menu_items_;

    /// the windows submenu (holds all windows added via addWindowToggle())
    QMenu* m_windows_;
  };

} //namespace

