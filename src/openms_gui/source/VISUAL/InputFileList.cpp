// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/InputFileList.h>
#include <ui_InputFileList.h>

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFileDialog.h>

#include <QtWidgets/QFileDialog>
#include <QApplication>
#include <QClipboard>
#include <QKeyEvent>
#include <QUrl>
#include <QMimeData>

#include <QString>
#include <QtWidgets/QFileDialog>
#include <QMessageBox>
#include <QListWidgetItem>

//#include <iostream>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    InputFileList::InputFileList(QWidget *parent) :
      QWidget(parent),
      ui_(new Ui::InputFileList)
    {
      ui_->setupUi(this);
      connect(ui_->add_button, SIGNAL(clicked()), this, SLOT(showFileDialog()));
      connect(ui_->edit_button, SIGNAL(clicked()), this, SLOT(editCurrentItem()));
      connect(ui_->remove_button, SIGNAL(clicked()), this, SLOT(removeSelected()));
      connect(ui_->remove_all_button, SIGNAL(clicked()), this, SLOT(removeAll()));
    }

    InputFileList::~InputFileList()
    {
      delete ui_;
    }

    void InputFileList::dragEnterEvent(QDragEnterEvent* e)
    {
      // file dropped from a window manager come as URLs
      if (e->mimeData()->hasUrls())
      {
        e->acceptProposedAction();
      }
    }

    void InputFileList::dropEvent(QDropEvent* e)
    {
      QStringList files;
      for (const QUrl& url : e->mimeData()->urls())
      {
        files  << url.toLocalFile();
      }
      addFiles_(files);
    }

    void InputFileList::dragMoveEvent(QDragMoveEvent* p_event)
    {
      // TODO allow filtering?
      //if (!p_event->mimeData()->hasFormat(MY_MIMETYPE))
      //{
      //  p_event->ignore();
      //  return;
      //}
      p_event->accept();
    }

    void InputFileList::keyPressEvent(QKeyEvent* e) {
      // when Ctrl-C is pressed, copy all selected files to clipboard as text
      if (e->matches(QKeySequence::Copy))
      {
        QStringList strings;
        QList<QListWidgetItem*> selected_items = ui_->input_file_list->selectedItems();
        foreach(QListWidgetItem * item, selected_items)
        {
          strings << item->text();
        }
        QApplication::clipboard()->setText(strings.join("\n"));
        e->accept(); // do not propagate upstream
      }
      // exit on escape (without saving the list)
      else if (e->key() == Qt::Key_Escape)
      {
        this->close();
      }
      // delete currently selected items
      else if (e->key() == Qt::Key_Delete)
      {
        removeSelected();
      }
    }

    void InputFileList::showFileDialog()
    {
      QStringList file_names = QFileDialog::getOpenFileNames(this, tr("Select input file(s)"), cwd_);
      addFiles_(file_names);
    }

    void InputFileList::removeSelected()
    {
      QList<QListWidgetItem*> selected_items = ui_->input_file_list->selectedItems();
      for (QListWidgetItem * item : selected_items)
      {
        ui_->input_file_list->takeItem(ui_->input_file_list->row(item));
      }
      updateCWD_();
    }

    void InputFileList::removeAll()
    {
      ui_->input_file_list->clear();
      updateCWD_();
    }

    void InputFileList::getFilenames(QStringList& files) const
    {
      files.clear();
      for (int i = 0; i < ui_->input_file_list->count(); ++i)
      {
        files.push_back(ui_->input_file_list->item(i)->text());
      }
    }

    StringList InputFileList::getFilenames() const
    {
      int nr_files = ui_->input_file_list->count();
      StringList files;
      for (int i = 0; i < nr_files; ++i)
      {
        files.push_back(ui_->input_file_list->item(i)->text());
      }
      return files;
    }

    void OpenMS::Internal::InputFileList::setFilenames(const QStringList& files)
    {
      addFiles_(files);
    }

    const QString& InputFileList::getCWD() const
    {
      return cwd_;
    }

    void OpenMS::Internal::InputFileList::setCWD(const QString& cwd, bool force)
    {
      if (force || (cwd_.isEmpty() && !cwd.isEmpty())) // do not set cwd_ as empty (does not help the user in browsing for files)
      {
        cwd_ = cwd;
      }
      emit updatedCWD(cwd_);
    }

    void InputFileList::editCurrentItem()
    {
      QListWidgetItem* item = ui_->input_file_list->currentItem();
      if (!item)
      {
        if (ui_->input_file_list->count() == 0)
        {
          return;
        }
        // use the first item if none is selected
        ui_->input_file_list->setCurrentItem(ui_->input_file_list->item(0));
        item = ui_->input_file_list->currentItem();
      }

      TOPPASInputFileDialog tifd(item->text());
      if (tifd.exec())
      {
        item->setText(tifd.getFilename());
        updateCWD_();
      }
    }

    void InputFileList::addFiles_(const QStringList& files)
    {
      if (!files.isEmpty())
      {
        ui_->input_file_list->addItems(files);
        setCWD(File::path(files.back()).toQString()); // emit the signal
      }
    }

    void OpenMS::Internal::InputFileList::updateCWD_()
    {
      QListWidgetItem* item = ui_->input_file_list->currentItem();
      // also update with empty, to ensure emitting the updatedCWD() signal
      setCWD(item ? item->text() : "", false);
    }

  } //namespace Internal
} //namspace OpenMS

