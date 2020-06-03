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
        return;
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

