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

#include <OpenMS/VISUAL/RecentFilesMenu.h>

#include <OpenMS/DATASTRUCTURES/String.h>
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
    if (!action) return;
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
