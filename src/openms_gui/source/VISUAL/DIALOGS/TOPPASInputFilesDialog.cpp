// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFilesDialog.h>
#include <ui_TOPPASInputFilesDialog.h>


#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFileDialog.h>

#include <OpenMS/SYSTEM/File.h>

#include <QtWidgets/QFileDialog>
#include <QApplication>
#include <QClipboard>
#include <QKeyEvent>
#include <QUrl>
#include <QMimeData>

#include <iostream>

namespace OpenMS
{
  TOPPASInputFilesDialog::TOPPASInputFilesDialog(const QStringList & list, const QString& cwd)
    : cwd_(cwd),
      ui_(new Ui::TOPPASInputFilesDialogTemplate)
  {
    ui_->setupUi(this);

    //input_file_list->setSortingEnabled(true);
    ui_->input_file_list->addItems(list);

    connect(ui_->ok_button, SIGNAL(clicked()), this, SLOT(accept()));
    connect(ui_->cancel_button, SIGNAL(clicked()), this, SLOT(reject()));
    connect(ui_->add_button, SIGNAL(clicked()), this, SLOT(showFileDialog()));
    connect(ui_->remove_button, SIGNAL(clicked()), this, SLOT(removeSelected()));
    connect(ui_->remove_all_button, SIGNAL(clicked()), this, SLOT(removeAll()));
    connect(ui_->edit_button, SIGNAL(clicked()), this, SLOT(editCurrentItem()));
    connect(ui_->up_button, SIGNAL(clicked()), this, SLOT(moveCurrentItem()));
    connect(ui_->down_button, SIGNAL(clicked()), this, SLOT(moveCurrentItem()));

    // allow dragging of filenames from OS window manager (Finder, Explorer etc)
    setAcceptDrops(true);
  }

  TOPPASInputFilesDialog::~TOPPASInputFilesDialog()
  {
    delete ui_;
  }

  void TOPPASInputFilesDialog::dragEnterEvent(QDragEnterEvent* e)
  {
    // file dropped from a window manager come as URLs
    if (e->mimeData()->hasUrls())
    {
      e->acceptProposedAction();
    }
  }

  void TOPPASInputFilesDialog::dropEvent(QDropEvent* e)
  {
    foreach (const QUrl& url, e->mimeData()->urls())
    {
      ui_->input_file_list->addItem(url.toLocalFile());
    }
  }

  void TOPPASInputFilesDialog::keyPressEvent(QKeyEvent* e) {
    // when Ctrl-C is pressed, copy all selected files to clipboard as text
    if (e->matches(QKeySequence::Copy))
    {
      QStringList strings;
      QList<QListWidgetItem*> selected_items = ui_->input_file_list->selectedItems();
      foreach (QListWidgetItem* item, selected_items)
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
  }

  void TOPPASInputFilesDialog::showFileDialog()
  {
    QStringList file_names = QFileDialog::getOpenFileNames(this,
                                                           tr("Select input file(s)"), 
                                                           cwd_);
    if (!file_names.isEmpty())
    {
      ui_->input_file_list->addItems(file_names);
      cwd_ = File::path(file_names.back()).toQString();
    }
  }

  void TOPPASInputFilesDialog::removeSelected()
  {
    QList<QListWidgetItem *> selected_items = ui_->input_file_list->selectedItems();
    foreach(QListWidgetItem * item, selected_items)
    {
      ui_->input_file_list->takeItem(ui_->input_file_list->row(item));
    }
  }

  void TOPPASInputFilesDialog::removeAll()
  {
    ui_->input_file_list->clear();
  }

  void TOPPASInputFilesDialog::getFilenames(QStringList& files) const
  {
    files.clear();
    for (int i = 0; i < ui_->input_file_list->count(); ++i)
    {
      files.push_back(ui_->input_file_list->item(i)->text());
    }
    if (ui_->flag_sort_list->isChecked())
      files.sort();
  }

  const QString& TOPPASInputFilesDialog::getCWD() const
  {
    return cwd_;
  }

  void TOPPASInputFilesDialog::editCurrentItem()
  {
    QListWidgetItem * item = ui_->input_file_list->currentItem();
    if (!item)
    {
      return;
    }
    TOPPASInputFileDialog tifd(item->text());
    if (tifd.exec())
    {
      item->setText(tifd.getFilename());
    }
  }

  void TOPPASInputFilesDialog::moveCurrentItem()
  {
    if (ui_->input_file_list->count() < 2)
    {
      return;
    }
    int row = ui_->input_file_list->currentRow();
    if (row < 0)
    {
      return;
    }

    if (QObject::sender() == ui_->up_button)     // move upwards
    {
      if (row == 0)
      {
        return;
      }
      QListWidgetItem * item = ui_->input_file_list->takeItem(row);
      ui_->input_file_list->insertItem(row - 1, item);
      ui_->input_file_list->setCurrentItem(item);
    }
    else if (QObject::sender() == ui_->down_button) // move downwards
    {
      if (row == ui_->input_file_list->count() - 1)
      {
        return;
      }
      QListWidgetItem * item = ui_->input_file_list->takeItem(row);
      ui_->input_file_list->insertItem(row + 1, item);
      ui_->input_file_list->setCurrentItem(item);
    }
  }

} // namespace
