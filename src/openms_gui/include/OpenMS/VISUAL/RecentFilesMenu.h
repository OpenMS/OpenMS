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

#include <OpenMS/KERNEL/StandardTypes.h>

#include <QMenu>
#include <QStringList>

#include <vector>

class QAction;

namespace OpenMS
{
  class Param;
  class String;

  /**
    @brief Manages recent files opened by the user and provides a QMenu to go with it


  */
  class RecentFilesMenu
    : public QObject
  {
    Q_OBJECT
  
  signals:
    /// when a recent file action item from the getMenu() was clicked
    void recentFileClicked(const String& filename);
    
  public:
    /// C'tor
    RecentFilesMenu(int max_entries = 15);

    /// sets a list of recent files (up to max_entries many -- see C'tor)
    void set(const QStringList& initial);

    /** 
        @brief Extracts all values from all elements in the param object and tries to interpret them as filenames
        If they exist, they will be used in the list of recent files. The name of the param items is ignored.

        @param filenames A Param object of which all values will be tested for being a filename
        @return The number of items which were successfully interpreted as filenames

    */
    unsigned setFromParam(const Param& filenames);

    /**
       @brief Convert current file list to Param. Their names are just numbers, starting at "0". The values are the filenames.
       @return Param object with name:value pairs
    */
    Param getAsParam() const;

    /// get a menu-pointer to an internal member which always contains the up-to-date recent items
    QMenu* getMenu();
    
    /// current list of recent files (most recent first)
    const QStringList& get() const;

  public slots:
    /// put a new recent file at the top (removing any duplicates in other positions); will update the QMenu
    void add(const String& filename);

  private slots:
    /// invoked by the QAction when it was clicked; emits recentFileClicked(String filename)
    void itemClicked_();

  private:
    /// updates the menu by synching text and and visibility of actions using the current list of recent files
    void sync_();

    /// holds the menu and the filenames (as QActions)
    QMenu recent_menu_;
    /// maximum of entries; adding more will delete the oldest one
    int max_entries_;
    /// list of the recently opened files actions (menu entries)
    QStringList recent_files_;
    /// .. and the actions to go with it
    std::vector<QAction*> recent_actions_;
  };

} //namespace

