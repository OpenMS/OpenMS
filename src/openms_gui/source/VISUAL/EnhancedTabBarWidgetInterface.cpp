// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/EnhancedTabBarWidgetInterface.h>

#include <OpenMS/VISUAL/EnhancedTabBar.h>

#include <QObject>

namespace OpenMS
{

  EnhancedTabBarWidgetInterface::EnhancedTabBarWidgetInterface()
  {
    /// every new window gets a new ID automatically
    static Int window_counter_ = getFirstWindowID();
    window_id_ = ++window_counter_;
  }

  EnhancedTabBarWidgetInterface::~EnhancedTabBarWidgetInterface()
  { // we cannot emit signals (since we cannot derive from QObject), so we let our member do it
    sp_.emitAboutToBeDestroyed(window_id_);
  }

  void EnhancedTabBarWidgetInterface::addToTabBar(EnhancedTabBar* const parent, const String& caption, const bool make_active_tab)
  {
    // use signal/slot to communicate, since directly storing the parent pointer for later access is dangerous (it may already be destroyed during program exit)
    QObject::connect(&this->sp_, &SignalProvider::aboutToBeDestroyed, parent, &EnhancedTabBar::removeId);
    parent->addTab(caption.toQString(), window_id_);
    if (make_active_tab)
    {
      parent->show(window_id_);
    }
  }

  Int EnhancedTabBarWidgetInterface::getWindowId() const
  {
    return window_id_;
  }

  /*static*/ Int EnhancedTabBarWidgetInterface::getFirstWindowID()
  {
    return 1234;
  }

}

