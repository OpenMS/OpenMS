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

#include <OpenMS/KERNEL/StandardTypes.h>

#include <QMenu>
#include <QStringList>

#include <vector>

class QAction;

namespace OpenMS
{
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

