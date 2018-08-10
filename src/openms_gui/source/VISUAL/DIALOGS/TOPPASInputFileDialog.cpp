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
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFileDialog.h>
#include <ui_TOPPASInputFileDialog.h>

#include <QtWidgets/QMessageBox>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QCompleter>
#include <QtWidgets/QDirModel>
#include <QtCore/QFileInfo>

#include <iostream>

namespace OpenMS
{
  TOPPASInputFileDialog::TOPPASInputFileDialog(const QString & file_name)
    : ui_(new Ui::TOPPASInputFileDialogTemplate)
  {
    ui_->setupUi(this);

    ui_->line_edit->setText(file_name);
    // disable completer for windows, causes crashes
#ifndef OPENMS_WINDOWSPLATFORM
    QCompleter * completer = new QCompleter(this);
    completer->setModel(new QDirModel(completer));
    ui_->line_edit->setCompleter(completer);
#endif
    connect(ui_->browse_button, SIGNAL(clicked()), this, SLOT(showFileDialog()));
    connect(ui_->ok_button, SIGNAL(clicked()), this, SLOT(checkValidity_()));
    connect(ui_->cancel_button, SIGNAL(clicked()), this, SLOT(reject()));
  }

  TOPPASInputFileDialog::~TOPPASInputFileDialog()
  {
    delete ui_;
  }

  void TOPPASInputFileDialog::showFileDialog()
  {
    QString file_name = QFileDialog::getOpenFileName(this, tr("Specify input file"), tr(""), tr(/*valid formats*/ ""));
    if (file_name != "")
    {
      ui_->line_edit->setText(file_name);
    }
  }

  QString TOPPASInputFileDialog::getFilename()
  {
    return ui_->line_edit->text();
  }

  void TOPPASInputFileDialog::checkValidity_()
  {
    // we ALLOW non-existing filenames (e.g. for FASTA files, which
    // are searched in other paths via OpenMS.ini:id_db_dir
    if (!fileNameValid(ui_->line_edit->text()))
    {
      QMessageBox::warning(nullptr, "Invalid file name", "Warning: filename does not exist!");
    }

    accept();
  }

  bool TOPPASInputFileDialog::fileNameValid(const QString & file_name)
  {
    QFileInfo fi(file_name);
    return fi.exists() && fi.isReadable() && (!fi.isDir());
  }

} // namespace
