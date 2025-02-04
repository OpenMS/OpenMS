// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/RecentFilesMenu.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/SYSTEM/File.h>

#include <QAction>

/*
#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>

#include <QAction>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
*/

using namespace std;

namespace OpenMS
{
  RecentFilesMenu::RecentFilesMenu(int max_entries)
    : recent_menu_("&Recent files"),
    max_entries_(max_entries),
    recent_files_()
  {
    // add hidden actions
    recent_actions_.resize(max_entries_);
    for (int i = 0; i < max_entries_; ++i)
    {
      recent_actions_[i] = recent_menu_.addAction("", this, &RecentFilesMenu::itemClicked_);
      recent_actions_[i]->setVisible(false);
    }
  }

  void RecentFilesMenu::set(const QStringList& initial)
  {
    recent_files_ = initial;
    recent_files_.removeDuplicates();
    while (recent_files_.size() > max_entries_)
    {
      recent_files_.removeLast();
    }
    sync_();
  }

  unsigned RecentFilesMenu::setFromParam(const Param& filenames)
  {
    QStringList rfiles;
    unsigned count{ 0 };
    for (Param::ParamIterator it = filenames.begin(); it != filenames.end(); ++it)
    {
      QString filename = String(it->value.toString()).toQString();
      if (File::exists(filename))
      {
        rfiles.append(filename);
        ++count;
      }
    }
    set(rfiles);
    return count;
  }

  Param RecentFilesMenu::getAsParam() const
  {
    Param p;
    int i{ 0 };
    for (const auto& f : recent_files_)
    {
      p.setValue(String(i), f.toStdString());
      ++i;
    }
    return p;
  }

  QMenu* RecentFilesMenu::getMenu()
  {
    return &recent_menu_;
  }

  const QStringList& RecentFilesMenu::get() const
  {
    return recent_files_;
  }

  void RecentFilesMenu::add(const String& filename)
  {
    // find out absolute path
    String tmp = File::absolutePath(filename);

    // remove the new file if already in the recent list and prepend it
    recent_files_.removeAll(tmp.toQString());
    recent_files_.prepend(tmp.toQString());

    // remove those files exceeding the defined number
    while (recent_files_.size() > max_entries_)
    {
      recent_files_.removeLast();
    }
    sync_();
  }

  void RecentFilesMenu::itemClicked_()
  {
    QAction* action = qobject_cast<QAction*>(sender());
    if (!action)
    {
      return;
    }
    String filename = String(action->text());
    emit recentFileClicked(filename);
  }

  void RecentFilesMenu::sync_()
  {
    for (int i = 0; i < max_entries_; ++i)
    {
      if (i < recent_files_.size())
      {
        recent_actions_[i]->setText(recent_files_[i]);
        recent_actions_[i]->setVisible(true);
      }
      else
      {
        recent_actions_[i]->setVisible(false);
      }
    }
  }

} //Namespace
