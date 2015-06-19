// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFilesDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFileDialog.h>

#include <QtGui/QFileDialog>
#include <QApplication>
#include <QClipboard>
#include <QKeyEvent>
#include <QUrl>

#include <iostream>

namespace OpenMS
{
  TOPPASInputFilesDialog::TOPPASInputFilesDialog(const QStringList & list)
  {
    setupUi(this);

    //input_file_list->setSortingEnabled(true);
    input_file_list->addItems(list);

    connect(ok_button, SIGNAL(clicked()), this, SLOT(accept()));
    connect(cancel_button, SIGNAL(clicked()), this, SLOT(reject()));
    connect(add_button, SIGNAL(clicked()), this, SLOT(showFileDialog()));
    connect(remove_button, SIGNAL(clicked()), this, SLOT(removeSelected()));
    connect(remove_all_button, SIGNAL(clicked()), this, SLOT(removeAll()));
    connect(edit_button, SIGNAL(clicked()), this, SLOT(editCurrentItem()));
    connect(up_button, SIGNAL(clicked()), this, SLOT(moveCurrentItem()));
    connect(down_button, SIGNAL(clicked()), this, SLOT(moveCurrentItem()));

    // allow dragging of filenames from OS window manager (Finder, Explorer etc)
    setAcceptDrops(true);
  }

  void TOPPASInputFilesDialog::dragEnterEvent(QDragEnterEvent* e)
  {
    // file dropped from a window manager come as URLs
    if (e->mimeData()->hasUrls()) {
      e->acceptProposedAction();
    }
  }

  void TOPPASInputFilesDialog::dropEvent(QDropEvent* e)
  {
    foreach (const QUrl& url, e->mimeData()->urls())
    {
      input_file_list->addItem(url.toLocalFile());
    }
  }

  void TOPPASInputFilesDialog::keyPressEvent(QKeyEvent* e) {
    // when Ctrl-C is pressed, copy all selected files to clipboard as text
    if (e->matches(QKeySequence::Copy))
    {
      QStringList strings;
      QList<QListWidgetItem*> selected_items = input_file_list->selectedItems();
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
    QStringList file_names = QFileDialog::getOpenFileNames(this, tr("Select input file(s)"), tr(""), tr(/*valid filetypes*/ ""));
    if (!file_names.isEmpty())
    {
      input_file_list->addItems(file_names);
    }
  }

  void TOPPASInputFilesDialog::removeSelected()
  {
    QList<QListWidgetItem *> selected_items = input_file_list->selectedItems();
    foreach(QListWidgetItem * item, selected_items)
    {
      input_file_list->takeItem(input_file_list->row(item));
    }
  }

  void TOPPASInputFilesDialog::removeAll()
  {
    input_file_list->clear();
  }

  void TOPPASInputFilesDialog::getFilenames(QStringList& files) const
  {
    files.clear();
    for (int i = 0; i < input_file_list->count(); ++i)
    {
      files.push_back(input_file_list->item(i)->text());
    }
    if (flag_sort_list->isChecked())
      files.sort();
  }

  void TOPPASInputFilesDialog::editCurrentItem()
  {
    QListWidgetItem * item = input_file_list->currentItem();
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
    if (input_file_list->count() < 2)
    {
      return;
    }
    int row = input_file_list->currentRow();
    if (row < 0)
    {
      return;
    }

    if (QObject::sender() == up_button)     // move upwards
    {
      if (row == 0)
      {
        return;
      }
      QListWidgetItem * item = input_file_list->takeItem(row);
      input_file_list->insertItem(row - 1, item);
      input_file_list->setCurrentItem(item);
    }
    else if (QObject::sender() == down_button) // move downwards
    {
      if (row == input_file_list->count() - 1)
      {
        return;
      }
      QListWidgetItem * item = input_file_list->takeItem(row);
      input_file_list->insertItem(row + 1, item);
      input_file_list->setCurrentItem(item);
    }
  }

} // namespace
