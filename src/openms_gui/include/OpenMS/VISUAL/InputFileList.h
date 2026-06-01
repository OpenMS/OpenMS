// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QWidget>

namespace Ui
{
  class InputFileList;
}

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief A widget shows a list of input files (i.e. existing files on a mounted drive),
             which allows adding/removing files and supports drag'n'drop from the window manager.

    */
    class InputFileList : public QWidget
    {
        Q_OBJECT

    public:
        /// C'tor
        explicit InputFileList(QWidget* parent = nullptr);
        ~InputFileList();
        /// support drag'n'drop of files from OS window manager
        void dragEnterEvent(QDragEnterEvent* e) override;
        /// support drag'n'drop of files from OS window manager
        void dropEvent(QDropEvent* e) override;
        void dragMoveEvent(QDragMoveEvent* pEvent) override;
        
        /// Stores the list of all filenames in the list widget in @p files
        void getFilenames(QStringList& files) const;
        /// Stores the list of all filenames in the list widget in @p files
        StringList getFilenames() const;
        /// Set the list of all filenames in the list widget
        void setFilenames(const QStringList& files);

        /// get the CWD (according to most recently added file)
        const QString& getCWD() const;
        /// set the current working directory (for opening files), but only if the current input list is not already populated. Use @p force to set the CWD in any case.
        void setCWD(const QString& cwd, bool force = false);

        /// support Ctrl+C to copy currently selected items to clipboard
        void keyPressEvent(QKeyEvent* e) override;

    public slots:

      /// Lets the user select files via a file dialog
      void showFileDialog();
      /// Removes all currently selected files from the list
      void removeSelected();
      /// Removes all files from the list
      void removeAll();
      /// Shows a TOPPASInputFileDialog which edits the current item
      void editCurrentItem();

    signals:
      /// emitted when a new file is added (by drag'n'drop or 'Add..' button)
      void updatedCWD(QString new_cwd);

    protected:
      /// add files to the list, and update 'cwd_' by using the path of the last filename
      void addFiles_(const QStringList& files);

      /// updates the CWD, based on the last file in the current list
      void updateCWD_();


      QString cwd_; ///< current working dir, i.e. the last position a file was added from

    private:
      Ui::InputFileList *ui_;
    };
  } // ns Internal
} // ns OpenMS

// this is required to allow parent widgets (auto UIC'd from .ui) to have a InputFileList member
using InputFileList = OpenMS::Internal::InputFileList;
