// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QObject>

namespace OpenMS
{
  class EnhancedTabBar;

  /**
    @brief provides a signal mechanism (by deriving from QObject) for classes which are not allowed to have signals themselves.

    This might be useful for EnhancedTabBarWidgetInterface, since that cannot derive from QObject due to the diamond star inheritance problem via its parent classes (e.g. PlotWidget).
    
    Diamond star problem:
      
      PlotWidget
       /       \
    ETBWI    QWidget
       -!      /
        QObject

     Thus, ETBWI cannot derive from QObject and needs to delegate its signaling duties to a SignalProvider.      

     Wrap all signals that are required in a function call and call these functions instead of emitting the signal directly.
     Connect the signal to a slot by using QObject::connect() externally somewhere.

  */
  class OPENMS_GUI_DLLAPI SignalProvider
    : public QObject
  {
    Q_OBJECT
  public:
    void emitAboutToBeDestroyed(int id)
    {
      emit aboutToBeDestroyed(id);
    }
  signals:
    void aboutToBeDestroyed(int id);
  };
  
  /**
    @brief Widgets that are placed into an EnhancedTabBar must implement this interface

    @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI EnhancedTabBarWidgetInterface
  {
  public:
    /// C'tor; creates a new ID;
    EnhancedTabBarWidgetInterface();
    /// Destructor (emits SignalProvider::aboutToBeDestroyed)
    virtual ~EnhancedTabBarWidgetInterface();

    /// adds itself to this tabbar and upon destruction removes itself again.
    /// Make sure the tabbar still exists when you call this function and this object is destroyed
    void addToTabBar(EnhancedTabBar* const parent, const String& caption, const bool make_active_tab = true);

    /// get the EnhancedTabBar unique window id
    Int getWindowId() const;

    /// the first object to be created will get this ID
    static Int getFirstWindowID();

  private:
    Int window_id_ { -1 };
    SignalProvider sp_; ///< emits the signal that the EnhancedTabBarWidgetInterface is about to be destroyed
  };
}  // namespace OpenMS

